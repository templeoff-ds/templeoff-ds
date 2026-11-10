# Построчный разбор parameter_study.py

Дополнительные проверки каждого настраиваемого числа; базовая модель не подменяется альтернативами.

[Исходный файл](../parameter_study.py) · [Общий маршрут чтения](../code-walkthrough.md)

Каждая физическая строка приведена отдельно. Если несколько действий соединены `;`, они выполняются слева направо. Длинная строка может переноситься экраном, но номер относится к исходному файлу.

## Интерфейсы функций

- **slope_result**, строка 22: регрессия ln d по t в window и R²; при нехватке точек возвращает None.
- **study**, строка 32: самостоятельная серия проверок меньших и больших настроек в constants-study.json.
- **save**, строка 42: локальный помощник сохранения одной фигуры или таблицы; использует путь из внешней функции.

## Строка за строкой

| № | Исходная строка | Что и зачем делает |
| --- | --- | --- |
| 1 | <code>&quot;&quot;&quot;Дополнительные эксперименты для защиты выбора численных констант.</code> | Документирующая строка модуля/функции. Объясняет назначение кода, сама не интегрирует и не меняет параметры. |
| 2 | — | Пустая строка: визуально разделяет блоки, при выполнении ничего не делает. |
| 3 | <code>Запуск: python parameter_study.py. Исходная модель и её результаты не меняются.</code> | Документирующая строка модуля/функции. Объясняет назначение кода, сама не интегрирует и не меняет параметры. |
| 4 | <code>Диапазоны альтернатив явно заданы в constants-study.json, а не спрятаны в коде.</code> | Документирующая строка модуля/функции. Объясняет назначение кода, сама не интегрирует и не меняет параметры. |
| 5 | <code>&quot;&quot;&quot;</code> | Документирующая строка модуля/функции. Объясняет назначение кода, сама не интегрирует и не меняет параметры. |
| 6 | <code>import ast</code> | Импортирует ast, но в текущем parameter_study.py это имя далее не используется. Он не влияет на численные эксперименты; AST реально применяется в отдельном build_constant_index.py. |
| 7 | <code>import json</code> | Импорт: сериализация и чтение JSON. Имена становятся доступны далее в этом модуле. |
| 8 | <code>import platform</code> | Импорт: фактическая версия Python. Имена становятся доступны далее в этом модуле. |
| 9 | <code>from pathlib import Path</code> | Импорт: работа с путями как объектами Path. Имена становятся доступны далее в этом модуле. |
| 10 | <code>from time import perf_counter</code> | Импорт: высокоточный монотонный таймер perf_counter. Имена становятся доступны далее в этом модуле. |
| 11 | <code>import numpy as np</code> | Импорт: векторные массивы, нормы и линейная алгебра (сокращение np). Имена становятся доступны далее в этом модуле. |
| 12 | <code>import scipy</code> | Импорт: версия библиотеки и научные вычисления. Имена становятся доступны далее в этом модуле. |
| 13 | <code>import matplotlib</code> | Импорт: настройка backend рисования. Имена становятся доступны далее в этом модуле. |
| 14 | <code>matplotlib.use(&#x27;Agg&#x27;)</code> | До pyplot включает серверный backend Agg; дополнительные рисунки не требуют GUI. |
| 15 | <code>import matplotlib.pyplot as plt</code> | Импорт: фигуры и оси под сокращением plt. Имена становятся доступны далее в этом модуле. |
| 16 | <code>from lorenz import Parameters, rhs, jacobian, equilibria</code> | Импорт: параметры, правая часть, Якобиан и аналитические равновесия. Имена становятся доступны далее в этом модуле. |
| 17 | <code>from numerical_methods import integrate, time_grid</code> | Импорт: собственные формулы шага, сетка и фиксированное интегрирование. Имена становятся доступны далее в этом модуле. |
| 18 | <code>from experiments import load_config, adaptive, write_csv</code> | Импорт: конфигурация, эксперименты и обёртки измерений/записи. Имена становятся доступны далее в этом модуле. |
| 19 | — | Пустая строка: визуально разделяет блоки, при выполнении ничего не делает. |
| 20 | <code>ROOT=Path(__file__).resolve().parent</code> | Находит папку проекта по расположению этого скрипта. |
| 21 | — | Пустая строка: визуально разделяет блоки, при выполнении ничего не делает. |
| 22 | <code>def slope_result(t,d,window):</code> | Определяет slope_result(t, d, window): регрессия ln d по t в window и R²; при нехватке точек возвращает None. Тело выполняется при вызове функции; отступы ниже задают его границы. |
| 23 | <code>    mask=(t&gt;=window[0])&amp;(t&lt;=window[1])&amp;(d&gt;0)</code> | Выбирает точки внутри window, у которых расстояние положительно. Для каждого t получается булево значение. |
| 24 | <code>    if np.count_nonzero(mask)&lt;2:</code> | Считает True в маске. При менее двух точках прямую с двумя коэффициентами оценить нельзя. |
| 25 | <code>        return dict(slope=None,r_squared=None,d_min=float(d.min()),d_max=float(d.max()))</code> | Возвращает отсутствующие наклон/R² и диапазон всего d. В этой ветви d_min/d_max относятся ко всему массиву, а не к пустому окну; это важно при интерпретации отказа. |
| 26 | <code>    a,b=np.polyfit(t[mask],np.log(d[mask]),1)</code> | На выбранных точках строит прямую по ln d; a — наклон, b — свободный член. |
| 27 | <code>    residual=np.sum((np.log(d[mask])-(a*t[mask]+b))**2)</code> | Считает SSE — сумму квадратов остатков прямой в логарифмическом пространстве. |
| 28 | <code>    total=np.sum((np.log(d[mask])-np.log(d[mask]).mean())**2)</code> | Считает SST — сумму квадратов отклонений ln d от своего среднего. |
| 29 | <code>    return dict(slope=float(a),r_squared=float(1-residual/total) if total&gt;0 else None,</code> | Возвращает наклон и R²; при нулевом SST показатель отсутствует, деления на ноль нет. |
| 30 | <code>                d_min=float(d[mask].min()),d_max=float(d[mask].max()))</code> | Добавляет минимальное и максимальное d только в допустимом окне (успешная ветвь). |
| 31 | — | Пустая строка: визуально разделяет блоки, при выполнении ничего не делает. |
| 32 | <code>def study():</code> | Определяет study(): самостоятельная серия проверок меньших и больших настроек в constants-study.json. Тело выполняется при вызове функции; отступы ниже задают его границы. |
| 33 | <code>    c=load_config();q=json.loads((ROOT/&#x27;constants-study.json&#x27;).read_text())</code> | Читает основную проверенную конфигурацию c и отдельный список альтернатив q; исходные параметры не перезаписываются. |
| 34 | <code>    out=ROOT/&#x27;results&#x27;/&#x27;constants_study&#x27;;out.mkdir(parents=True,exist_ok=True)</code> | Создаёт отдельную папку constants_study, чтобы не смешивать дополнительные CSV с исходной серией. |
| 35 | <code>    p=Parameters(c[&#x27;sigma&#x27;],c[&#x27;rho&#x27;],c[&#x27;beta&#x27;]);u=np.array(c[&#x27;initial_state&#x27;]);fun=lambda t,x:rhs(t,x,p)</code> | Фиксирует классические параметры/начальный вектор и функцию скорости с этими параметрами. |
| 36 | <code>    t=time_grid(c[&#x27;t_start&#x27;],c[&#x27;t_end&#x27;],c[&#x27;step&#x27;])</code> | Строит общую основную сетку. |
| 37 | <code>    base=integrate(fun,u,c[&#x27;t_start&#x27;],c[&#x27;t_end&#x27;],c[&#x27;step&#x27;])</code> | Считает базовый RK4 с параметрами по умолчанию; явный method не указан, integrate по умолчанию использует rk4. |
| 38 | <code>    pert=integrate(fun,u+np.array(c[&#x27;perturbation&#x27;]),c[&#x27;t_start&#x27;],c[&#x27;t_end&#x27;],c[&#x27;step&#x27;])</code> | Считает вторую RK4-траекторию с исходным возмущением. |
| 39 | <code>    d=np.linalg.norm(pert.y-base.y,axis=1)</code> | Получает расстояния пары по всем моментам основной сетки. |
| 40 | <code>    ref=adaptive(fun,u,c[&#x27;t_start&#x27;],c[&#x27;t_end&#x27;],c[&#x27;reference_rtol&#x27;],c[&#x27;reference_atol&#x27;])</code> | Считает строгий RK45 как численное сравнение, а не аналитическое решение. |
| 41 | <code>    result={&#x27;versions&#x27;:{&#x27;python&#x27;:platform.python_version(),&#x27;numpy&#x27;:np.__version__,&#x27;scipy&#x27;:scipy.__version__},&#x27;baseline&#x27;:c,&#x27;study_config&#x27;:q}</code> | Создаёт результат с версиями, исходной конфигурацией и полным набором альтернатив для воспроизводимости. |
| 42 | <code>    def save(name,rows):</code> | Определяет save(name, rows): локальный помощник сохранения одной фигуры или таблицы; использует путь из внешней функции. Тело выполняется при вызове функции; отступы ниже задают его границы. |
| 43 | <code>        result[name]=rows;write_csv(out/(name+&#x27;.csv&#x27;),rows)</code> | Сохраняет одну таблицу и в общем словаре, и отдельным CSV с переданным именем. |
| 44 | <code>    rows=[]</code> | Создаёт новый пустой список для следующей таблицы; не дополняет случайно результаты предыдущей подсерии. |
| 45 | <code>    for key,values in q[&#x27;model_variants&#x27;].items():</code> | Выбирает очередной изменяемый параметр модели и список его значений. |
| 46 | <code>        for value in values:</code> | Перебирает все альтернативы именно этого параметра. |
| 47 | <code>            pp=Parameters(**{**dict(sigma=p.sigma,rho=p.rho,beta=p.beta),key:value})</code> | Создаёт новый словарь классических параметров и заменяет только key; ** распаковывает словарь в именованные аргументы Parameters. Старый p не меняется. |
| 48 | <code>            ss=adaptive(lambda t,x:rhs(t,x,pp),u,c[&#x27;t_start&#x27;],c[&#x27;t_end&#x27;],c[&#x27;reference_rtol&#x27;],c[&#x27;reference_atol&#x27;])</code> | Решает изменённую модель строгим RK45; lambda передаёт новый pp, а не базовый p. |
| 49 | <code>            ys=ss.sol(t).T;tail=ys[t&gt;=c[&#x27;transient_end&#x27;]]</code> | Получает состояния на общей сетке и хвост после transient_end. Изменяются физические параметры, поэтому это другая задача, не ошибка относительно прежней траектории. |
| 50 | <code>            rows.append(dict(parameter=key,value=value,max_state_norm=float(np.max(np.linalg.norm(ys,axis=1))),</code> | Начинает строку с именем параметра, значением и максимальной длиной состояния за интервал. |
| 51 | <code>                             tail_std_x=float(np.std(tail[:,0])),</code> | Сохраняет стандартное отклонение x на хвосте как описательную характеристику, не критерий хаоса. |
| 52 | <code>                             equilibria_count=len(equilibria(pp)),max_real_eigen_Cplus=float(np.max(np.real(np.linalg.eigvals(jacobian(equilibria(pp)[&#x27;C+&#x27;],pp))))) if pp.rho&gt;1 else None))</code> | Сохраняет число равновесий и максимальную Re(lambda) в C+, только если C+ существует; иначе None. |
| 53 | <code>    save(&#x27;model_parameters&#x27;,rows)</code> | Передаёт заполненную таблицу помощнику save: она попадёт в result["model_parameters"] и файл model_parameters.csv. |
| 54 | <code>    rows=[]</code> | Создаёт новый пустой список для следующей таблицы; не дополняет случайно результаты предыдущей подсерии. |
| 55 | <code>    for initial in q[&#x27;initial_variants&#x27;]:</code> | Перебирает альтернативные начальные векторы при неизменной модели. |
| 56 | <code>        sol=adaptive(fun,initial,c[&#x27;t_start&#x27;],c[&#x27;t_end&#x27;],c[&#x27;reference_rtol&#x27;],c[&#x27;reference_atol&#x27;])</code> | Для каждого нового старта решает соответствующую задачу Коши. |
| 57 | <code>        ys=sol.sol(t).T</code> | Оценивает это решение на общей временной сетке, форма (M,3). |
| 58 | <code>        rows.append(dict(initial=str(initial),initial_speed=float(np.linalg.norm(fun(t[0],initial))),max_norm=float(np.max(np.linalg.norm(ys,axis=1))),final_norm=float(np.linalg.norm(ys[-1]))))</code> | Записывает старт, модуль исходной скорости, максимальную и конечную длину состояния. Старт O даёт нулевую скорость. |
| 59 | <code>    save(&#x27;initial_states&#x27;,rows)</code> | Передаёт заполненную таблицу помощнику save: она попадёт в result["initial_states"] и файл initial_states.csv. |
| 60 | <code>    # Общий численный эталон; для каждого h сравнение идёт на собственной сетке sol.t.</code> | Комментарий автора: Общий численный эталон; для каждого h сравнение идёт на собственной сетке sol.t.. Комментарий не выполняется. |
| 61 | <code>    rows=[]</code> | Создаёт новый пустой список для следующей таблицы; не дополняет случайно результаты предыдущей подсерии. |
| 62 | <code>    for end in q[&#x27;short_horizons&#x27;]:</code> | Выбирает очередную длину короткого интервала. |
| 63 | <code>        errors=[]</code> | Создаёт отдельный список пар шаг–ошибка для данного горизонта. |
| 64 | <code>        for h in sorted(c[&#x27;steps&#x27;],reverse=True):</code> | Перебирает исходные h от больших к меньшим. |
| 65 | <code>            sol=integrate(fun,u,c[&#x27;t_start&#x27;],end,h,&#x27;rk4&#x27;) if np.isclose(round((end-c[&#x27;t_start&#x27;])/h)*h,end-c[&#x27;t_start&#x27;]) else None</code> | Считает RK4 только если этот интервал согласуется с h. Иначе устанавливает None; в данной дополнительной серии скрытый остаточный шаг не добавляется. |
| 66 | <code>            if sol is None:continue</code> | Пропускает несогласованную пару интервал–шаг; она не участвует в наклоне как нулевая ошибка. |
| 67 | <code>            error=float(np.max(np.linalg.norm(sol.y-ref.sol(sol.t).T,axis=1)))</code> | Вычисляет максимальное расхождение RK4 с RK45 на узлах этого h и горизонта. |
| 68 | <code>            errors.append((h,error))</code> | Добавляет пару (h,E) для регрессии. |
| 69 | <code>        chosen=errors[-3:]</code> | Берёт три мелкие успешные сетки этого горизонта. |
| 70 | <code>        order=float(np.polyfit(np.log([r[0] for r in chosen]),np.log([r[1] for r in chosen]),1)[0])</code> | Строит log–log наклон порядка по выбранным парам. |
| 71 | <code>        rows.append(dict(end=end,order_last_three=order,error_h001=next(e for h,e in errors if h==c[&#x27;step&#x27;])))</code> | Записывает горизонт, порядок и ошибку при основном h. next требует, чтобы основной шаг присутствовал; набор альтернатив рассчитан на исходный config. |
| 72 | <code>    save(&#x27;short_intervals&#x27;,rows)</code> | Передаёт заполненную таблицу помощнику save: она попадёт в result["short_intervals"] и файл short_intervals.csv. |
| 73 | <code>    horizons=[]</code> | Создаёт список характеристик разных полных горизонтов. |
| 74 | <code>    for end in q[&#x27;horizons&#x27;]:</code> | Перебирает заданные времена окончания, включая горизонт длиннее основного. |
| 75 | <code>        if end &lt;= c[&#x27;t_end&#x27;]:</code> | Для горизонтов внутри основной сетки можно использовать уже рассчитанный массив. |
| 76 | <code>            error=float(np.max(np.linalg.norm(base.y[t&lt;=end]-ref.sol(t[t&lt;=end]).T,axis=1)))</code> | Берёт соответствующий начальный участок базового RK4 и сравнивает с эталоном на тех же моментах. |
| 77 | <code>            distance=float(d[np.searchsorted(t,end)])</code> | searchsorted находит первый индекс t≥end; для текущих согласованных горизонтов это нужный узел. Произвольный несеточный end потребовал бы отдельного уточнения метрики. |
| 78 | <code>        else:</code> | Альтернативная ветвь непосредственно предшествующего if: выполняется, когда его условие ложно. |
| 79 | <code>            longer=integrate(fun,u,c[&#x27;t_start&#x27;],end,c[&#x27;step&#x27;])</code> | Для более длинного горизонта самостоятельно продолжает расчёт от исходного u0 на всём новом интервале. |
| 80 | <code>            longer_p=integrate(fun,u+np.array(c[&#x27;perturbation&#x27;]),c[&#x27;t_start&#x27;],end,c[&#x27;step&#x27;])</code> | Аналогично решает возмущённую задачу на удлинённом интервале. |
| 81 | <code>            longer_ref=adaptive(fun,u,c[&#x27;t_start&#x27;],end,c[&#x27;reference_rtol&#x27;],c[&#x27;reference_atol&#x27;])</code> | Отдельно получает численный эталон для длинного интервала. |
| 82 | <code>            error=float(np.max(np.linalg.norm(longer.y-longer_ref.sol(longer.t).T,axis=1)))</code> | Вычисляет максимальное расхождение с новым эталоном. Позднее значение не объявляется строгой точной ошибкой. |
| 83 | <code>            distance=float(np.linalg.norm(longer.y[-1]-longer_p.y[-1]))</code> | Считает конечное расстояние удлинённой пары близких траекторий. |
| 84 | <code>        n=int(round((end-c[&#x27;t_start&#x27;])/c[&#x27;step&#x27;]))</code> | Выводит количество фиксированных шагов из длины интервала и h. |
| 85 | <code>        horizons.append(dict(end=end,max_rk4_error=error,distance_at_end=distance,n_steps=n,nfev=4*n))</code> | Записывает время окончания, ошибки/расстояние, число шагов и точное для RK4 4N вызовов. |
| 86 | <code>    save(&#x27;horizons&#x27;,horizons)</code> | Сохраняет таблицу horizons отдельно и в общем result, включая стоимость длинных расчётов. |
| 87 | <code>    save(&#x27;distance_thresholds&#x27;,[dict(threshold=threshold,first_crossing=float(t[np.flatnonzero(d&gt;=threshold)[0]]) if np.any(d&gt;=threshold) else None) for threshold in q[&#x27;distance_thresholds&#x27;]])</code> | Для каждого диагностического порога находит первое достижение d≥порог, если оно есть. Не пересчитывает динамику: меняется только интерпретация одной кривой d(t). |
| 88 | <code>    short_grid=time_grid(c[&#x27;t_start&#x27;],c[&#x27;short_end&#x27;],min(c[&#x27;steps&#x27;]))</code> | Создаёт общую мелкую выходную сетку для всех пар допусков. |
| 89 | <code>    control=adaptive(fun,u,c[&#x27;t_start&#x27;],c[&#x27;short_end&#x27;],q[&#x27;control_rtol&#x27;],q[&#x27;control_atol&#x27;],&#x27;DOP853&#x27;).sol(short_grid).T</code> | Считает DOP853 с контрольными допусками и оценивает его на этой сетке. |
| 90 | <code>    rows=[]</code> | Создаёт новый пустой список для следующей таблицы; не дополняет случайно результаты предыдущей подсерии. |
| 91 | <code>    for rtol,atol in q[&#x27;tolerance_pairs&#x27;]:</code> | Берёт по очереди пары rtol/atol, которые меняются совместно. |
| 92 | <code>        seconds=[]</code> | Создаёт список длительностей для этой пары допусков. |
| 93 | <code>        for _ in range(c[&#x27;timing_repeats&#x27;]):</code> | Повторяет измерение столько раз, сколько указано в основной конфигурации. |
| 94 | <code>            begin=perf_counter();sol=adaptive(fun,u,c[&#x27;t_start&#x27;],c[&#x27;short_end&#x27;],rtol,atol)</code> | Засекает время и решает короткую задачу RK45 с данной парой допусков. |
| 95 | <code>            ys=sol.sol(short_grid).T;seconds.append(perf_counter()-begin)</code> | Оценивает решение на общей сетке; длительность включает интерполяцию, а не только вызов solve_ivp. |
| 96 | <code>        rows.append(dict(rtol=rtol,atol=atol,max_error=float(np.max(np.linalg.norm(ys-control,axis=1))),nfev=sol.nfev,n_steps=len(sol.t)-1,seconds_median=float(np.median(seconds))))</code> | Сохраняет допуски, расхождение с контролем, реальные внутренние счётчики и медиану времени. |
| 97 | <code>    save(&#x27;tolerances&#x27;,rows)</code> | Передаёт заполненную таблицу помощнику save: она попадёт в result["tolerances"] и файл tolerances.csv. |
| 98 | <code>    rows=[]</code> | Создаёт новый пустой список для следующей таблицы; не дополняет случайно результаты предыдущей подсерии. |
| 99 | <code>    for delta in q[&#x27;perturbations&#x27;]:</code> | Перебирает величины начального сдвига первой координаты. |
| 100 | <code>        perturbed_u=u.copy();perturbed_u[0]+=delta</code> | Обязательно копирует u, затем прибавляет delta только к первому элементу копии. Без copy испортился бы общий базовый вектор последующих экспериментов. |
| 101 | <code>        ys=integrate(fun,perturbed_u,c[&#x27;t_start&#x27;],c[&#x27;t_end&#x27;],c[&#x27;step&#x27;]).y</code> | Интегрирует RK4 для конкретного возмущения. |
| 102 | <code>        dist=np.linalg.norm(ys-base.y,axis=1);reached=np.flatnonzero(dist&gt;=1)</code> | Получает расстояния с base и индексы достижения порога 1. |
| 103 | <code>        rows.append(dict(requested_delta=delta,actual_delta=float(perturbed_u[0]-u[0]),first_above_one=float(t[reached[0]]) if reached.size else None,**slope_result(t,dist,c[&#x27;slope_window&#x27;])))</code> | Сохраняет запрошенную и фактически представленную delta, первый момент достижения и статистику наклона. ** распаковывает возвращённый slope_result словарь в строку. |
| 104 | <code>    save(&#x27;perturbations&#x27;,rows)</code> | Передаёт заполненную таблицу помощнику save: она попадёт в result["perturbations"] и файл perturbations.csv. |
| 105 | <code>    save(&#x27;windows&#x27;,[dict(start=a,end=b,**slope_result(t,d,[a,b])) for a,b in q[&#x27;slope_windows&#x27;]])</code> | Перебирает окна на одной и той же кривой d, не выбирая новый старт/метод. У каждого окна своя регрессия. |
| 106 | <code>    save(&#x27;discard_times&#x27;,[dict(discard=begin,retained_fraction=float(np.mean(t&gt;=begin)),retained_points=int(np.count_nonzero(t&gt;=begin)),std_x=float(np.std(base.y[t&gt;=begin,0]))) for begin in q[&#x27;discard_times&#x27;]])</code> | Сравнивает удаление разных начальных отрезков: mean от булевой маски — доля сохранённых узлов, count_nonzero — их число, std — разброс x на хвосте. |
| 107 | <code>    # Эффект числа последних точек, а не подбор лучшего наклона.</code> | Комментарий автора: Эффект числа последних точек, а не подбор лучшего наклона.. Комментарий не выполняется. |
| 108 | <code>    existing=json.loads((ROOT/&#x27;results&#x27;/&#x27;summary.json&#x27;).read_text())</code> | Читает исходный summary.json для анализа выбора количества регрессионных точек. |
| 109 | <code>    if existing[&#x27;config&#x27;] != c:</code> | Сверяет сохранённую конфигурацию с текущей, чтобы не сравнить новые настройки со старыми ошибками. |
| 110 | <code>        raise ValueError(&#x27;Конфигурация изменилась: сначала повторите python main.py&#x27;)</code> | Просит сначала пересчитать основной эксперимент, если конфигурации разошлись. Проверка выполняется здесь, после предыдущих подсерий; это не откат уже записанных дополнительных CSV. |
| 111 | <code>    rows=[]</code> | Создаёт новый пустой список для следующей таблицы; не дополняет случайно результаты предыдущей подсерии. |
| 112 | <code>    for method in [&#x27;euler&#x27;,&#x27;rk4&#x27;]:</code> | Рассматривает выбор точек отдельно для Эйлера и RK4. |
| 113 | <code>        good=[r for r in existing[&#x27;convergence&#x27;] if r[&#x27;method&#x27;]==method and r[&#x27;status&#x27;]==&#x27;ok&#x27;]</code> | Фильтрует только успешные исходные строки нужного метода. |
| 114 | <code>        for count in q[&#x27;fit_points&#x27;]:</code> | Пробует разные количества последних точек из q. |
| 115 | <code>            if count&gt;len(good):continue</code> | Если запросил больше точек, чем есть успешных, пропускает вариант. Например, у Эйлера грубый отказ уменьшает доступное число. |
| 116 | <code>            chosen=good[-count:]</code> | Выделяет count последних строк, то есть мелких h в исходном порядке. |
| 117 | <code>            rows.append(dict(method=method,point_count=count,h_max=max(r[&#x27;h&#x27;] for r in chosen),h_min=min(r[&#x27;h&#x27;] for r in chosen),order=float(np.polyfit(np.log([r[&#x27;h&#x27;] for r in chosen]),np.log([r[&#x27;max_error&#x27;] for r in chosen]),1)[0])))</code> | Сохраняет количество, диапазон h и наклон. Таблица показывает зависимость оценки от поднабора, не подставляет заранее нужный p. |
| 118 | <code>    save(&#x27;fit_points&#x27;,rows)</code> | Передаёт заполненную таблицу помощнику save: она попадёт в result["fit_points"] и файл fit_points.csv. |
| 119 | <code>    rows=[]</code> | Создаёт новый пустой список для следующей таблицы; не дополняет случайно результаты предыдущей подсерии. |
| 120 | <code>    for repeats in q[&#x27;timing_repeat_counts&#x27;]:</code> | Пробует несколько чисел повторов замера времени. |
| 121 | <code>        durations=[]</code> | Обнуляет список длительностей для нового размера серии. |
| 122 | <code>        for _ in range(repeats):</code> | Запускает ровно выбранное число измерений. |
| 123 | <code>            begin=perf_counter();integrate(fun,u,c[&#x27;t_start&#x27;],c[&#x27;t_end&#x27;],c[&#x27;step&#x27;]);durations.append(perf_counter()-begin)</code> | Измеряет только собственное RK4-интегрирование; создаваемые здесь траектории не используются для подмены основных результатов. |
| 124 | <code>        rows.append(dict(repeats=repeats,median_seconds=float(np.median(durations)),min_seconds=float(min(durations)),max_seconds=float(max(durations))))</code> | Сохраняет медиану и наблюдаемый min/max времени. Это не формальный доверительный интервал. |
| 125 | <code>    save(&#x27;timings&#x27;,rows)</code> | Передаёт заполненную таблицу помощнику save: она попадёт в result["timings"] и файл timings.csv. |
| 126 | <code>    state=np.array([2.3,-1.4,8.1]);rows=[]</code> | Выбирает ту же нетривиальную точку, что и в тесте Якобиана, и начинает новую таблицу. |
| 127 | <code>    for eps in q[&#x27;jacobian_eps&#x27;]:</code> | Перебирает положительные приращения для центральной разности. |
| 128 | <code>        fd=np.column_stack([(fun(0,state+eps*e)-fun(0,state-eps*e))/(2*eps) for e in np.eye(3)])</code> | По каждому единичному направлению считает разностный вектор скорости, собирает три столбца численного J. |
| 129 | <code>        rows.append(dict(eps=eps,max_entry_error=float(np.max(np.abs(fd-jacobian(state,p))))))</code> | Вычитает аналитический J, берёт поэлементные абсолютные ошибки и наибольшую из всех девяти. Это не спектральная матричная норма. |
| 130 | <code>    save(&#x27;jacobian_difference&#x27;,rows)</code> | Передаёт заполненную таблицу помощнику save: она попадёт в result["jacobian_difference"] и файл jacobian_difference.csv. |
| 131 | <code>    rows=[]</code> | Создаёт новый пустой список для следующей таблицы; не дополняет случайно результаты предыдущей подсерии. |
| 132 | <code>    for threshold in q[&#x27;eigen_thresholds&#x27;]:</code> | Перебирает технические пороги различения нулевых Re(lambda). |
| 133 | <code>        for name,state in equilibria(p).items():</code> | Для каждого порога проверяет каждое равновесие той же модели. |
| 134 | <code>            es=np.real(np.linalg.eigvals(jacobian(state,p)))</code> | Получает вещественные части спектра J в данном равновесии. |
| 135 | <code>            label=&#x27;unstable&#x27; if np.any(es&gt;threshold) else &#x27;stable&#x27; if np.all(es &lt; -threshold) else &#x27;inconclusive&#x27;</code> | Повторяет критерий unstable/stable/inconclusive с текущим порогом, не меняя истинный спектр. |
| 136 | <code>            rows.append(dict(threshold=threshold,equilibrium=name,classification=label))</code> | Сохраняет тройку порог–равновесие–численная классификация. |
| 137 | <code>    save(&#x27;eigen_thresholds&#x27;,rows)</code> | Передаёт заполненную таблицу помощнику save: она попадёт в result["eigen_thresholds"] и файл eigen_thresholds.csv. |
| 138 | <code>    rows=[]</code> | Создаёт новый пустой список для следующей таблицы; не дополняет случайно результаты предыдущей подсерии. |
| 139 | <code>    for rtol,atol in q[&#x27;grid_tolerances&#x27;]:</code> | Перебирает пары допусков проверки делимости сетки, не допуски RK45. |
| 140 | <code>        # Первый случай имеет только двоичное округление, второй — намеренное несовпадение.</code> | Комментарий автора: Первый случай имеет только двоичное округление, второй — намеренное несовпадение.. Комментарий не выполняется. |
| 141 | <code>        for label,start,end,h in [(&#x27;binary_roundoff&#x27;,0,.3,.1),(&#x27;near_noninteger&#x27;,0,1+1e-12,.01),(&#x27;noninteger&#x27;,0,1+1e-7,.01)]:</code> | Задаёт три демонстрационных случая: обычное округление десятичных дробей, очень малое несовпадение длины, явно большее несовпадение. |
| 142 | <code>            n=round((end-start)/h)</code> | Округляет теоретическое число шагов для текущего случая. |
| 143 | <code>            rows.append(dict(rtol=rtol,atol=atol,case=label,length_disagreement=float(abs(n*h-(end-start))),accepted=bool(np.isclose(n*h,end-start,rtol=rtol,atol=atol))))</code> | Сохраняет фактическое несовпадение длины и решение isclose о его допустимости. accepted не означает, что ОДУ стало точнее. |
| 144 | <code>    save(&#x27;grid_tolerances&#x27;,rows)</code> | Передаёт заполненную таблицу помощнику save: она попадёт в result["grid_tolerances"] и файл grid_tolerances.csv. |
| 145 | <code>    # Ошибки, с которыми сравниваются численные пороги автоматических тестов.</code> | Комментарий автора: Ошибки, с которыми сравниваются численные пороги автоматических тестов.. Комментарий не выполняется. |
| 146 | <code>    rows=[]</code> | Создаёт новый пустой список для следующей таблицы; не дополняет случайно результаты предыдущей подсерии. |
| 147 | <code>    exp_error=abs(integrate(lambda t,x:x,[1],0,1,.01).y[-1,0]-np.e)</code> | Измеряет ошибку RK4 на точном тесте exp(t); это независимая проверка порога соответствующего unit-теста. |
| 148 | <code>    nonauto=abs(integrate(lambda t,x:t*x,[1],0,1,.01).y[-1,0]-np.exp(.5))</code> | Измеряет ошибку неавтономного теста с точным exp(t²/2). |
| 149 | <code>    check=integrate(fun,u,0,2,.001);control_dense=adaptive(fun,u,0,2,q[&#x27;control_rtol&#x27;],q[&#x27;control_atol&#x27;],&#x27;DOP853&#x27;)</code> | Считает мелкий RK4 и контроль DOP853 для короткого Лоренца. |
| 150 | <code>    short_error=float(np.max(np.linalg.norm(check.y-control_dense.sol(check.t).T,axis=1)))</code> | Определяет фактический максимум их расхождения на сетке RK4. |
| 151 | <code>    for name,error,threshold in [(&#x27;rk4_exp&#x27;,exp_error,3e-10),(&#x27;rk4_nonautonomous&#x27;,nonauto,1e-8),(&#x27;rk4_lorenz&#x27;,short_error,2e-6)]:</code> | Связывает три измеренные ошибки с фактически установленными в тестах порогами. |
| 152 | <code>        rows.append(dict(test=name,measured_error=float(error),threshold=threshold,threshold_over_error=float(threshold/error)))</code> | Сохраняет отношение порога к измеренной ошибке: значение больше 1 оставляет запас, но не доказывает оптимальность порога. |
| 153 | <code>    save(&#x27;test_error_margins&#x27;,rows)</code> | Передаёт заполненную таблицу помощнику save: она попадёт в result["test_error_margins"] и файл test_error_margins.csv. |
| 154 | <code>    plt.rcParams.update({&#x27;font.family&#x27;:&#x27;DejaVu Sans&#x27;,&#x27;font.size&#x27;:10})</code> | Настраивает кириллический шрифт для дополнительных изображений. |
| 155 | <code>    rows=[]</code> | Создаёт новый пустой список для следующей таблицы; не дополняет случайно результаты предыдущей подсерии. |
| 156 | <code>    for dpi in q[&#x27;dpi_variants&#x27;]:</code> | Перебирает разрешения, включая намеренный пример ниже требования. |
| 157 | <code>        fig,ax=plt.subplots(figsize=(6.4,4.8));ax.plot(t,base.y[:,0]);ax.set(xlabel=&#x27;t&#x27;,ylabel=&#x27;x&#x27;,title=f&#x27;Контроль вывода: dpi={dpi}&#x27;)</code> | Создаёт одинаковый временной график x(t) для сравнения dpi; меняется подпись разрешения, не численная траектория. |
| 158 | <code>        begin=perf_counter();filename=out/f&#x27;dpi_{dpi}.png&#x27;;fig.savefig(filename,dpi=dpi);elapsed=perf_counter()-begin;plt.close(fig)</code> | Измеряет сохранение PNG, получает его имя, закрывает фигуру. elapsed относится к записи/отрисовке, а не интегрированию. |
| 159 | <code>        from PIL import Image</code> | Импорт: чтение сохранённого изображения через Pillow. Имена становятся доступны далее в этом модуле. |
| 160 | <code>        with Image.open(filename) as im:width,height=im.size</code> | Открывает сохранённую картинку и читает реальную ширину/высоту в пикселях. |
| 161 | <code>        rows.append(dict(dpi=dpi,width_px=width,height_px=height,bytes=filename.stat().st_size,save_seconds=elapsed))</code> | Сохраняет dpi, размеры, вес файла и время записи. |
| 162 | <code>    save(&#x27;image_resolution&#x27;,rows)</code> | Передаёт заполненную таблицу помощнику save: она попадёт в result["image_resolution"] и файл image_resolution.csv. |
| 163 | <code>    rows=[]</code> | Создаёт новый пустой список для следующей таблицы; не дополняет случайно результаты предыдущей подсерии. |
| 164 | <code>    for nframes in q[&#x27;frame_variants&#x27;]:</code> | Перебирает возможные количества кадров, не создавая девять новых GIF. |
| 165 | <code>        idx=np.unique(np.linspace(1,len(t),nframes).astype(int))-1</code> | Выбирает уникальные целые индексы кадров; −1 переводит номера концов 1..M в индексы элементов 0..M−1. |
| 166 | <code>        for fps in q[&#x27;fps_variants&#x27;]:</code> | Для данного набора кадров пробует разные скорости воспроизведения. |
| 167 | <code>            rows.append(dict(frames=len(idx),fps=fps,duration_seconds=float(len(idx)/fps),max_simulation_time_gap=float(np.max(np.diff(t[idx]))),all_integrator_points_still_in_line=True))</code> | Вычисляет номинальную длительность и максимальный разрыв модельного времени между маркерами. True подчёркивает, что линия содержит все накопленные шаги интегратора. |
| 168 | <code>    save(&#x27;animation_sampling&#x27;,rows)</code> | Передаёт заполненную таблицу помощнику save: она попадёт в result["animation_sampling"] и файл animation_sampling.csv. |
| 169 | <code>    # Сводные рисунки не затрагивают исходные графики.</code> | Комментарий автора: Сводные рисунки не затрагивают исходные графики.. Комментарий не выполняется. |
| 170 | <code>    fig,axes=plt.subplots(2,2,figsize=(12,9))</code> | Создаёт четыре панели сводного рисунка компромиссов. |
| 171 | <code>    r=result[&#x27;tolerances&#x27;];axes[0,0].loglog([x[&#x27;nfev&#x27;] for x in r],[x[&#x27;max_error&#x27;] for x in r],&#x27;o-&#x27;);axes[0,0].set(xlabel=&#x27;Вызовы правой части&#x27;,ylabel=&#x27;Ошибка на [0,2]&#x27;,title=&#x27;Допуски: точность и стоимость&#x27;)</code> | Рисует точность как функцию стоимости nfev в log–log шкале для пар допусков. |
| 172 | <code>    r=result[&#x27;perturbations&#x27;];axes[0,1].semilogx([x[&#x27;actual_delta&#x27;] for x in r if x[&#x27;actual_delta&#x27;]&gt;0],[x[&#x27;first_above_one&#x27;] if x[&#x27;first_above_one&#x27;] is not None else np.nan for x in r if x[&#x27;actual_delta&#x27;]&gt;0],&#x27;o-&#x27;);axes[0,1].set(xlabel=&#x27;Фактическое возмущение&#x27;,ylabel=&#x27;Первый t с d≥1&#x27;,title=&#x27;Меньшая δ откладывает расхождение&#x27;)</code> | Показывает фактическую delta и первое время достижения d≥1; непредставимое нулевое возмущение исключается из логарифмической оси, отсутствие достижения отображается пропуском NaN только в рисунке. |
| 173 | <code>    r=result[&#x27;windows&#x27;];axes[1,0].bar(range(len(r)),[x[&#x27;slope&#x27;] for x in r]);axes[1,0].set_xticks(range(len(r)),[f&quot;{x[&#x27;start&#x27;]}–{x[&#x27;end&#x27;]}&quot; for x in r]);axes[1,0].set(xlabel=&#x27;Окно регрессии&#x27;,ylabel=&#x27;Наклон&#x27;,title=&#x27;Зависимость оценки от окна&#x27;)</code> | Столбцами сравнивает наклоны для разных окон, подписывает каждое окно. |
| 174 | <code>    r=result[&#x27;jacobian_difference&#x27;];axes[1,1].loglog([x[&#x27;eps&#x27;] for x in r],[max(x[&#x27;max_entry_error&#x27;],np.finfo(float).tiny) for x in r],&#x27;o-&#x27;);axes[1,1].set(xlabel=&#x27;Приращение ε&#x27;,ylabel=&#x27;max ошибки J&#x27;,title=&#x27;Разностный Якобиан: слишком малое ε вредно&#x27;)</code> | Показывает ошибку разностного J по epsilon. tiny — нижнее положительное нормальное число float, используется только для допустимого логарифма нулевого значения на графике. |
| 175 | <code>    for ax in axes.flat:ax.grid(alpha=.25)</code> | Включает слабую сетку на всех четырёх осях через плоский обход axes.flat. |
| 176 | <code>    fig.tight_layout();fig.savefig(out/&#x27;constant_tradeoffs.png&#x27;,dpi=c[&#x27;dpi&#x27;]);plt.close(fig)</code> | Выравнивает панели, сохраняет сводный PNG с основным dpi и закрывает фигуру. |
| 177 | <code>    (out/&#x27;study.json&#x27;).write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding=&#x27;utf-8&#x27;)</code> | Записывает все таблицы, настройки и версии в study.json. |
| 178 | <code>    print(&#x27;Результаты проверки констант:&#x27;,out)</code> | Печатает папку завершённого исследования. |
| 179 | <code>    return result</code> | Возвращает тот же результат для возможного программного использования. |
| 180 | — | Пустая строка: визуально разделяет блоки, при выполнении ничего не делает. |
| 181 | <code>if __name__==&#x27;__main__&#x27;:study()</code> | Выполняет study только при прямом запуске скрипта. |
