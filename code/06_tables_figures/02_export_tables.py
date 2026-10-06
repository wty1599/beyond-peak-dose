"""Render released display tables as Markdown without re-estimating any value."""
from pathlib import Path
import csv

def cell(value):
    return value.replace('|','\\|').replace('\r','').replace('\n','<br>')

def main():
    destination=Path('results/private/table_previews')
    destination.mkdir(parents=True,exist_ok=True)
    for source in sorted(Path('results/aggregate').glob('Table_*_display.csv')):
        with source.open(encoding='utf-8-sig',newline='') as f:
            rows=list(csv.reader(f))
        if not rows:continue
        lines=['| '+' | '.join(map(cell,rows[0]))+' |',
               '| '+' | '.join(['---']*len(rows[0]))+' |']
        lines.extend('| '+' | '.join(map(cell,row))+' |' for row in rows[1:])
        (destination/(source.stem+'.md')).write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print('Released table previews written locally; no patient data were read.')

if __name__=='__main__':main()
