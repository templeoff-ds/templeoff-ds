"""Автоматическая вставка реальных результатов в русскоязычный отчёт."""
from pathlib import Path
import os
import json

def number(x):
    return '—' if x=='' or x is None else f'{x:.7g}'

def build_report(summary,output):
    root=Path(__file__).resolve().parent;out=Path(output).resolve();c=summary['config']
    text=['### 7.1. Настройки и среда\n', '```json\n'+json.dumps(c,ensure_ascii=False,indent=2)+'\n```\n',
          'Версии: '+', '.join(k+' '+v for k,v in summary['versions'].items())+'.\n',
          '### 7.2. Метрики на общей сетке\n',
          '| Метод | max ошибки [0,T] | Ошибка в T | max на коротком интервале | Время, с (медиана) | Шаги | Вызовы f | Точки вывода |',
          '| --- | --- | --- | --- | --- | --- | --- | --- |']
    for m in summary['metrics']:
        text.append('| '+' | '.join([m['method'],*[number(m[k]) for k in ['max_error','final_error','short_max_error','seconds_median']],str(m['n_steps']),str(m['nfev']),str(m['n_output'])])+' |')
    ref=summary['reference']
    text.extend(['\nДолгосрочная ошибка в этой таблице — расхождение с численным эталоном, а не достоверная оценка ошибки относительно точной траектории. Ранжирование методов делаем по короткому интервалу.\n',
                 f"Эталон RK45: rtol={ref['rtol']:g}, atol={ref['atol']:g}; {ref['n_steps']} принятых шагов, {ref['nfev']} вызовов f. Максимальное расхождение эталона с DOP853 на короткой контрольной сетке: **{number(ref['short_dop853_disagreement'])}**; с более строгим RK45: **{number(ref['short_tighter_rk45_disagreement'])}**. Это согласование, не строгая граница точности.\n",
                 '### 7.3. Исходная серия шагов\n','| Метод | h | Максимальная ошибка | Конечная ошибка | Порядок к предыдущему шагу | Статус |','| --- | --- | --- | --- | --- | --- |'])
    for r in summary['convergence']:
        text.append('| '+' | '.join([r['method'],number(r['h']),number(r['max_error']),number(r['final_error']),number(r['observed_order']),r['status']])+' |')
    text.append('\nnonfinite означает переполнение: конечное решение для этой пары «метод–шаг» не получено, поэтому ошибка и порядок не присваиваются. Это ожидаемый исследовательский результат, а не успешная проверка конечности этого прогона.\n')
    text.append('Наклон регрессии по трём самым мелким успешным шагам **исходной** серии: '+', '.join(k+': '+number(v) for k,v in summary['orders_last_three'].items())+'.\n')
    text.extend(['### 7.4. Дополнительная серия Эйлера (исходные результаты не заменены)\n','| h | Максимальная ошибка |','| --- | --- |'])
    for r in summary['euler_refinement']:text.append(f"| {number(r['h'])} | {number(r['max_error'])} |")
    text.append(f"\nНаклон по трём самым мелким дополнительным шагам: **{number(summary['euler_refinement_order'])}**. Приближение к единице демонстрирует асимптотический первый порядок. Дополнительные шаги заданы явно в конфигурации.\n")
    ss=summary['sensitivity']
    text.extend(['### 7.5. Чувствительность\n',f"Фактическое начальное расстояние: {number(ss['initial_distance'])}; расстояние в конце: {number(ss['final_distance'])}. Первый узел, где d≥1: {number(ss['first_distance_above_one'])} (если «—», порог не достигнут).\n",
                 f"Регрессия ln d(t)=a·t+b на окне {ss['window']}: a={number(ss['slope'])}, b={number(ss['intercept'])}, R²={number(ss['r_squared'])}. При h/2 наклон a={number(ss['fine_slope'])}. **Оба числа — конечновременные оценки для выбранного возмущения, не точный максимальный показатель Ляпунова.**\n",
                 '### 7.6. Равновесия и собственные значения\n','| Точка | Координаты | Собственные значения | Локальная устойчивость |','| --- | --- | --- | --- |'])
    for eq in summary['equilibria']:
        vals=[]
        for v in eq['eigenvalues']:vals.append(number(v['real'])+(f" {v['imag']:+.7g}i" if v['imag'] else ''))
        text.append('| '+' | '.join([eq['name'],', '.join(map(number,eq['state'])),'; '.join(vals),eq['stability']])+' |')
    text.append('\nКонечность всех основных траекторий (Эйлер/RK4/RK45 на основном шаге, эталон и возмущённый RK4): **'+str(summary['all_standard_trajectories_finite'])+'**. Это не включает заведомо неудачный крупный шаг из отдельной серии.\n')
    text.append('### 7.7. Визуальные результаты\n')
    captions={'lorenz_attractor.png':'Трёхмерная траектория после переходного процесса','time_series.png':'Временные зависимости x,y,z','phase_portraits.png':'Три фазовые проекции','methods_comparison.png':'Сравнение методов и поточечные расхождения','step_dependence.png':'Ошибка в зависимости от шага','sensitivity.png':'Чувствительность и зависимость оценки от шага','equilibria.png':'Равновесия и спектры Якобиана','animation.gif':'Анимация построения двух траекторий'}
    template=(root/'report.template.md').read_text(encoding='utf-8')
    def make(base):
        figures=[]
        for filename,caption in captions.items():
            if (out/filename).exists():figures.append('!['+caption+']('+Path(os.path.relpath(out/filename,base)).as_posix()+')\n\n'+caption+'.\n')
        content=template.replace('{{NUMERICAL_RESULTS}}','\n\n'.join(text+figures))
        for doc in ['constants-explained.md','constants-results.md','constants-index.md','math-glossary.md','code-walkthrough.md']:
            content=content.replace(']('+doc+')',']('+Path(os.path.relpath(root/doc,base)).as_posix()+')')
        return content
    (out/'report.md').write_text(make(out),encoding='utf-8')
    if out==root/'results':(root/'report.md').write_text(make(root),encoding='utf-8')
