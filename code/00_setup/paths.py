"""Read the JSON-compatible YAML path configuration. No credentials are stored here."""
from pathlib import Path
import json
import os

def release_path(key):
    config = Path(os.environ.get('BEYOND_PATHS_CONFIG', 'config/paths.yml'))
    if not config.is_file():
        raise RuntimeError('Create config/paths.yml from paths.example.yml.')
    values = json.loads(config.read_text(encoding='utf-8'))
    value = values.get(key)
    if not isinstance(value, str) or not value or value.startswith('REPLACE_'):
        raise RuntimeError(f'Configure {key}.')
    return Path(value)

def private_output(stage):
    p = release_path('private_output') / stage
    p.mkdir(parents=True, exist_ok=True)
    return p
