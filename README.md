# Starterator
Starterator is a program used to compare and analyze start sites of similar genes that have been
grouped together by phamerator.

Originally written by Marissa Pacey under the direction of Dan Russell.
It is currently under the maintenance and care of Chris Shaffer.

The current stable version can be found on github at [SEA-PHAGES/starterator](https://github.com/SEA-PHAGES/starterator);
any recent development can be found on github at [cdshaffer/starterator](https://github.com/cdshaffer/starterator).

## Database Backends (CLI)

Starterator CLI supports both MySQL (PyMySQL) and SQLite.

- MySQL mode (default): uses `DB_HOST`, `DB_USER`, `DB_PASSWORD`, `DB_NAME`, and `DB_PORT`.
- SQLite mode: set `DB_SQLITE_PATH` to an absolute path of a SQLite database file.

If `DB_SQLITE_PATH` is set and non-empty, SQLite mode is always used (it takes precedence over MySQL env vars).

### Local examples

```bash
# SQLite mode
DB_SQLITE_PATH=/absolute/path/actino_draft_v636.db \
python3 -m starterator.starterate --get-phams

# MySQL mode (default when DB_SQLITE_PATH is unset)
DB_HOST=127.0.0.1 DB_USER=root DB_PASSWORD=secret DB_NAME=actinophamerator \
python3 -m starterator.starterate -n 230 -j True
```

`--get-phams` returns only relevant phams with at least 2 genes that pass Starterator validation:
phage status is not `unknown`, genome sequence does not contain `N`, and each retained gene passes `has_valid_start()`.

### Docker example (SQLite mode)

```bash
docker run --rm \
  --user "$(id -u):$(id -g)" \
  -e DB_SQLITE_PATH=/db/actino_draft_v636.db \
  -e STARTERATOR_CONFIG_DIR=/app/config \
  -v "$PWD/starterator-out/reports:/app/reports:z" \
  -v "$PWD/starterator-out/alignments:/tmp/starterator_temp:z" \
  -v "$PWD/starterator-out/proteins:/app/starterator/Proteins:z" \
  -v "$PWD/starterator-out/config:/app/config:z" \
  -v "$HOME/Downloads:/db:ro,z" \
  starterator:latest --get-phams
```
