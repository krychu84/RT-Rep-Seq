#!/usr/bin/env bash
# RT-Rep-Seq installer: conda environment + configuration + `RT-Rep-Seq` command.
#
#   ./install.sh --db /path/to/metabuli/database          use an existing Metabuli database
#   ./install.sh --download-db gtdb226 --db-dir /data/db   download a pre-built database (large!)
#   ./install.sh                                           environment only; set the database later
#
# Options:
#   --env-name NAME      conda environment name (default: rt-rep-seq)
#   --db PATH            existing Metabuli database folder (contains db.parameters + taxonomy/)
#   --download-db NAME   gtdb226 (378 GB, recommended) | gtdb232 (744 GB) | refseq_standard (111 GB)
#                        | hrgm2 (85 GB) | hrom (42 GB)  -- Metabuli pre-built databases
#   --db-dir DIR         where to put a downloaded database (default: ~/databases/metabuli)
#   --bin-dir DIR        where to put the RT-Rep-Seq command (default: ~/.local/bin)
set -euo pipefail
export CONDA_ALWAYS_YES=true   # no interactive "Confirm changes" prompts

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="rt-rep-seq"; DB=""; DOWNLOAD=""; DB_DIR="$HOME/databases/metabuli"; BIN_DIR="$HOME/.local/bin"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --env-name) ENV_NAME="$2"; shift 2;;
    --db) DB="$2"; shift 2;;
    --download-db) DOWNLOAD="$2"; shift 2;;
    --db-dir) DB_DIR="$2"; shift 2;;
    --bin-dir) BIN_DIR="$2"; shift 2;;
    -h|--help) sed -n '2,16p' "$0"; exit 0;;
    *) echo "unknown option: $1 (see ./install.sh --help)"; exit 1;;
  esac
done

say()  { printf '\n\033[1;31m==>\033[0m \033[1m%s\033[0m\n' "$*"; }
ok()   { printf '    \033[32m%s\033[0m\n' "$*"; }
die()  { printf '\n\033[1;31mERROR:\033[0m %s\n' "$*" >&2; exit 1; }

say "1/5  Checking prerequisites"
[[ "$(uname -s)" == "Linux" ]] || echo "    note: RT-Rep-Seq is developed for Linux sequencing PCs (CPU/RAM caps use Linux features)"
YES=()
if command -v mamba >/dev/null 2>&1; then CONDA=mamba; YES=(-y)
elif command -v conda >/dev/null 2>&1; then CONDA=conda
else die "conda not found. Install Miniforge first: https://github.com/conda-forge/miniforge (then open a new terminal)"; fi
CONDA_BASE="$(conda info --base)"
ok "using $CONDA ($CONDA_BASE)"
command -v systemd-run >/dev/null 2>&1 && ok "systemd-run found (RAM cap can be enforced)" \
  || echo "    note: systemd-run not found -- the 50% RAM cap cannot be enforced (CPU cap still applies)"

say "2/5  Creating conda environment '$ENV_NAME' (Python, NumPy, SciPy, Biopython, Matplotlib, psutil, Metabuli)"
if [[ -d "$CONDA_BASE/envs/$ENV_NAME" && ! -d "$CONDA_BASE/envs/$ENV_NAME/conda-meta" ]]; then
  die "$CONDA_BASE/envs/$ENV_NAME exists but is not a complete conda environment (interrupted install?). Remove it: rm -rf $CONDA_BASE/envs/$ENV_NAME"
fi
if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  ok "environment exists -- updating"
  "$CONDA" env update "${YES[@]}" -n "$ENV_NAME" -f "$REPO/environment.yml" --prune
else
  "$CONDA" env create "${YES[@]}" -n "$ENV_NAME" -f "$REPO/environment.yml"
fi
PREFIX="$(conda env list | awk -v n="$ENV_NAME" '$1==n {print $NF}')"
PY="$PREFIX/bin/python"; METABULI="$PREFIX/bin/metabuli"
[[ -x "$PY" ]] || die "python not found in $PREFIX"
[[ -x "$METABULI" ]] || die "metabuli not found in $PREFIX"
"$PY" - <<EOF || die "Python packages could not be imported"
import sys; sys.path.insert(0, "$REPO")
import numpy, scipy, Bio, matplotlib, psutil, polyleven, torchlex
print("    packages OK: numpy", numpy.__version__, "| scipy", scipy.__version__, "| torchlex (bundled)", torchlex.__version__)
EOF
MB_VER="$(ls "$PREFIX"/conda-meta/metabuli-*.json 2>/dev/null | head -1 | xargs -r basename | cut -d- -f2)"
ok "metabuli ${MB_VER:-installed}: $METABULI"

say "3/5  Metabuli database"
if [[ -n "$DOWNLOAD" ]]; then
  case "$DOWNLOAD" in
    gtdb232) URL="https://steineggerlab.s3.amazonaws.com/metabuli/gtdb232.tar.gz"; SIZE=744;;
    gtdb226|refseq_standard|hrgm2|hrom) URL="https://opendata.mmseqs.org/metabuli/$DOWNLOAD.tar.gz"
      SIZE=$(case "$DOWNLOAD" in gtdb226) echo 378;; refseq_standard) echo 111;; hrgm2) echo 85;; hrom) echo 42;; esac);;
    *) die "unknown database '$DOWNLOAD' (gtdb226 | gtdb232 | refseq_standard | hrgm2 | hrom)";;
  esac
  mkdir -p "$DB_DIR"
  FREE=$(df -BG --output=avail "$DB_DIR" | tail -1 | tr -dc 0-9)
  (( FREE > SIZE * 2 )) || die "$DB_DIR has ${FREE} GB free; downloading + unpacking $DOWNLOAD needs about $((SIZE * 2)) GB"
  ok "downloading $DOWNLOAD (~${SIZE} GB) to $DB_DIR -- this can take hours"
  ( cd "$DB_DIR" && wget -c "$URL" && tar -xzf "$DOWNLOAD.tar.gz" && rm -f "$DOWNLOAD.tar.gz" )
  DB="$(dirname "$(find "$DB_DIR" -maxdepth 3 -name db.parameters -newer "$REPO/install.sh" | head -1)")"
fi
if [[ -n "$DB" ]]; then
  DB="$(cd "$DB" && pwd)"
  [[ -f "$DB/db.parameters" && -d "$DB/taxonomy" ]] || die "$DB does not look like a Metabuli database (db.parameters + taxonomy/ expected)"
  ok "database: $DB"
else
  echo "    no database given -- set it later with:  ./install.sh --db /path/to/database"
  echo "    (pre-built databases: https://jaebeom-kim.github.io/metabuli-doc/databases/new-database/ ; gtdb226 recommended)"
fi

say "4/5  Writing configuration (assets/rt-rep-seq.config.json)"
CFG="$REPO/assets/rt-rep-seq.config.json"
"$PY" - "$CFG" "$REPO/assets/rt-rep-seq.config.example.json" "$DB" "$METABULI" <<'EOF'
import json, os, sys
cfg_path, example, db, metabuli = sys.argv[1:5]
cfg = json.load(open(cfg_path)) if os.path.exists(cfg_path) else json.load(open(example))
if db:
    cfg['gtdb_db'] = db
cfg['metabuli_bin'] = metabuli
cfg['metabuli_env'] = ''
json.dump(cfg, open(cfg_path, 'w'), indent=2)
print('    gtdb_db      =', cfg.get('gtdb_db'))
print('    metabuli_bin =', metabuli)
print('    CPU / RAM budget = %d%% / %d%% of this computer' % (cfg.get('cpu_frac', 0.5) * 100, cfg.get('mem_frac', 0.5) * 100))
EOF

say "5/5  Installing the RT-Rep-Seq command"
mkdir -p "$BIN_DIR"
cat > "$BIN_DIR/RT-Rep-Seq" <<EOF
#!/usr/bin/env bash
# RT-Rep-Seq launcher (written by $REPO/install.sh)
export PATH="$PREFIX/bin:\$PATH"
exec "$PY" "$REPO/RT-Rep-Seq" "\$@"
EOF
chmod +x "$BIN_DIR/RT-Rep-Seq" "$REPO/RT-Rep-Seq"
ok "command: $BIN_DIR/RT-Rep-Seq"
case ":$PATH:" in
  *":$BIN_DIR:"*) ;;
  *) echo "export PATH=\"$BIN_DIR:\$PATH\"" >> "$HOME/.bashrc"
     echo "    added $BIN_DIR to PATH in ~/.bashrc -- open a new terminal (or: source ~/.bashrc)";;
esac

say "Done."
cat <<EOF
    Start an analysis:   cd /path/to/MinKNOW/<run>/<sample>/<flowcell_run>/fastq_pass
                         RT-Rep-Seq
    Live view again:     RT-Rep-Seq --status          Stop:  RT-Rep-Seq --stop
    Results:             ~/Desktop/RT-Rep-Seq_<run>_<date_time>/latest/report.html
EOF
[[ -n "$DB" ]] || printf '\n    \033[33mStill needed:\033[0m a Metabuli database -> ./install.sh --db /path/to/database\n'
