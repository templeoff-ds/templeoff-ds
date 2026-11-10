"""Точка входа: python main.py. Все пути разрешаются от каталога проекта."""
import argparse
from pathlib import Path
from experiments import load_config, run_experiments, DEFAULT_CONFIG
from visualization import render_figures, render_animation
from reporting import build_report

ROOT=Path(__file__).resolve().parent

def main():
    parser=argparse.ArgumentParser(description='Численное исследование системы Лоренца')
    parser.add_argument('--config',type=Path,default=DEFAULT_CONFIG)
    parser.add_argument('--output',type=Path,default=ROOT/'results')
    parser.add_argument('--skip-animation',action='store_true',help='Пропустить только пересоздание GIF')
    args=parser.parse_args();config=load_config(args.config)
    print('Численные эксперименты...',flush=True)
    data=run_experiments(config,args.output)
    print('Графики...',flush=True);render_figures(data,config,args.output)
    if not args.skip_animation:
        print('GIF-анимация...',flush=True);render_animation(data,config,args.output)
    build_report(data['summary'],args.output)
    print('Готово:',args.output.resolve())
    print('Наблюдаемые порядки (три малых шага):',data['summary']['orders_last_three'])

if __name__=='__main__':main()
