# -*- coding: utf-8 -*-
# DIM: na (генератор .tex з .md; ручки 0). К5.5 за нарізкою §K5.6.
"""
md_to_tex — compact_en.md / theory_en.md → .tex (ОДНОКОЛОНКОВО).

★ЧОМУ ГЕНЕРАТОР, А НЕ .tex РУКАМИ: написаний руками .tex — це ДРУГА ПРАВДА
про той самий текст. Через тиждень правка піде в .md, а PDF лишиться старим —
і ніхто не помітить, бо ніхто не звіряє .tex із .md рядок у рядок. Одне
джерело (.md) → .tex генерується. Правило-10: генероване руками не правити.

★ОДНОКОЛОНКОВО — слово автора (2026-07-17, дослівно): «двоколонкова не зручна
для швидкого читання, у нас не енциклопедія». Наслідок карбовано в нарізці:
компакт = ~3-3.7 стор., «один листочок» НЕ досягнутий, бо ~1100 із 1873 слів —
кваліфікатори, а §12 їх різати забороняє. Формат поступився правилу.

★★ГЕЙТ-ГЛІФ (несучий): будь-який не-ASCII символ, якого НЕМА в мапі, ВАЛИТЬ
генератор. ЧОМУ: невідомий гліф у LaTeX не кричить — він ТИХО зникає або стає
порожнім квадратом у PDF. Тобто «PDF зібрався» виглядало б однаково у світі
«все на місці» і у світі «з формули випав ≺». Символ, що зник із формули, —
це інша формула. ⟹ не знаємо символ — не випускаємо файл.

⚠★СТАТУС ВИХОДУ: .tex тут НЕ КОМПІЛЮЄТЬСЯ — у системі немає ні TeX, ні
pandoc (перевірено). ⟹ згенерований .tex — ЩЕ НЕ ПЕРЕВІРЕНИЙ артефакт:
відсутність тулчейна не робить його правильним. Поки він не зібраний у PDF
хоч десь (Overleaf/машина автора) — це КЛЕЙМ, не результат, і карбувати його
як готовий заборонено.

ЗАПУСК: cd active-v10.2 && python src/tools/md_to_tex.py [--check]
"""
import re
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

ROOT = Path(__file__).resolve().parents[2]
EXP = ROOT / "prime-metric-theory" / "export"

# ─────────────────────── МАПА ГЛІФІВ ──────────────────────────────────────────
# ★Кожен запис — рішення, не автоматика: символ, якого нема тут, валить гейт.
#
# ★★ІНДЕКСИ/СТЕПЕНІ — ПРОГОНАМИ, НЕ ПОСИМВОЛЬНО. Дві збірки поспіль довели,
#   що посимвольна мапа тут бреше:
#     «η⁻¹» → ^-^1     (інша формула, мовчки)
#     «K₀₂» → K_0_2    (! Double subscript — упало аж у pdflatex)
#   Причина одна: індекс «₀₂» — це ОДИН індекс із двох цифр, а не два індекси.
#   ⟹ знімаємо ПРОГІН символів у {…} до того, як до них дійде посимвольна мапа.
SUB = {"₀": "0", "₁": "1", "₂": "2", "₃": "3", "₋": "-"}
SUP = {"¹": "1", "²": "2", "³": "3", "⁴": "4", "⁻": "-", "ᵀ": r"\mathsf{T}"}
SUB_RUN = re.compile("[" + "".join(SUB) + "]+")
SUP_RUN = re.compile("[" + "".join(SUP) + "]+")


def _scripts(s):
    s = SUB_RUN.sub(lambda mo: "_{" + "".join(SUB[c] for c in mo.group(0)) + "}", s)
    return SUP_RUN.sub(lambda mo: "^{" + "".join(SUP[c] for c in mo.group(0)) + "}", s)

MATH = {
    "⟹": r"\Longrightarrow", "⟺": r"\Longleftrightarrow", "→": r"\to",
    "≺": r"\prec", "⊥": r"\perp", "∈": r"\in", "∧": r"\wedge",
    "∖": r"\setminus", "≠": r"\neq", "≤": r"\leq", "≥": r"\geq",
    "×": r"\times", "·": r"\cdot", "−": "-", "±": r"\pm",
    "Ω": r"\Omega", "η": r"\eta", "φ": r"\varphi", "ω": r"\omega",
    "𝟙": r"\mathbf{1}", 
     "ℝ": r"\mathbb{R}", "∎": r"\qed",
    "⊕": r"\oplus", "≡": r"\equiv", "∞": r"\infty", "√": r"\sqrt{}",
    "Δ": r"\Delta", "Φ": r"\Phi", "ε": r"\varepsilon", "π": r"\pi",
    "θ": r"\theta", "λ": r"\lambda", "κ": r"\kappa", "≅": r"\cong",
    "∥": r"\parallel", "∀": r"\forall", "∃": r"\exists", "⊂": r"\subset",
    "…": r"\dots", "α": r"\alpha", "β": r"\beta", "γ": r"\gamma",
    "∝": r"\propto", "≈": r"\approx", "≼": r"\preccurlyeq",
}
# ★∎ → \blacksquare, НЕ \qed: amsthm визначає \qed сам, і власний \newcommand
#   дає «Command \qed already defined» — спіймано ПЕРШОЮ ЖЕ збіркою.
MATH["∎"] = r"\blacksquare"


def _m(macro):
    """★Макрос, приклеєний до літери, — це ІНШИЙ макрос: «×p» дало «\\timesp»
    ⟹ Undefined control sequence ×7 на першій збірці. TeX завершує ім'я команди
    на першому не-літерному символі, тому після літерного макросу ОБОВ'ЯЗКОВИЙ
    роздільник. Ставимо пробіл (у математиці він нічого не займає).
    ★Цього НЕ бачив ні гейт-гліф (символи всі відомі), ні структурна звірка
    (оточення збалансовані) — лише компілятор."""
    return macro + " " if macro[:1] == "\\" and macro[-1:].isalpha() else macro
TEXT = {
    "§": r"\S{}", "★": r"$\bigstar$", "–": "--", "—": "---",
    "’": "'", "‘": "`", "“": "``", "”": "''", "…": r"\dots{}",
    "«": "``", "»": "''", "✓": r"$\checkmark$", "✗": r"$\times$",
    "⚠": r"$\bigstar$", "≈": r"$\approx$", "ä": r"\"{a}", "é": r"\'{e}",
    # ★† (U+2020) — ПОЗНАЧКА «акт записаний внутрішньо, публічної проби нема».
    #  Додано за виміром гейта-гліфа: без цього рядка знак ТИХО ЗНИКАВ БИ в PDF,
    #  і цитата S888† читалась би як звичайне посилання — тобто ремонт, зроблений
    #  саме проти тихої обіцянки, сам би зникнув по дорозі до читача. Гейт спіймав
    #  це ДО збірки й відмовився випускати: рівно те, заради чого він і стоїть.
    "†": r"$\dagger$",
}

# ★ОГОЛОШЕНИЙ ВИНЯТОК ГЕЙТА-ГЛІФ: КИРИЛИЦЯ.
#   theory_en.md несе цитати автора ДОСЛІВНО українською (+ підрядковий
#   переклад у дужках). Це НЕ недогляд перекладу: провенанс — це САМЕ СЛОВО
#   автора, і підмінити його англійським переказом означало б втратити те,
#   заради чого цитата стоїть. ⟹ .tex мусить УМІТИ кирилицю (T2A+babel),
#   а не викидати її.
#   ★Виняток ОГОЛОШЕНИЙ і ЛІЧЕНИЙ, не тихий: гейт друкує, скільки кириличних
#   символів пропущено. Мовчазний виняток тихо росте у відмичку (S904).
CYRILLIC = re.compile(r"[Ѐ-ӿ]")
# слова, що в математиці мусять стояти прямим шрифтом
MATH_WORDS = [("so", r"\mathfrak{so}"), ("Pf", r"\operatorname{Pf}"),
              ("det", r"\det"), ("ln", r"\ln"), ("dim", r"\dim"),
              ("span", r"\operatorname{span}"), ("Lie", r"\operatorname{Lie}"),
              ("diag", r"\operatorname{diag}"), ("Ad", r"\operatorname{Ad}"),
              ("iff", r"\text{ iff }")]


def _esc(s):
    """LaTeX-екранування ТЕКСТУ (не математики)."""
    for a, b in (("\\", r"\textbackslash{}"), ("&", r"\&"), ("%", r"\%"),
                 ("$", r"\$"), ("#", r"\#"), ("_", r"\_"), ("{", r"\{"),
                 ("}", r"\}"), ("~", r"\textasciitilde{}"), ("^", r"\textasciicircum{}")):
        s = s.replace(a, b)
    return s


def _math(s):
    """Вміст backtick'ів → математика."""
    for a, b in (("\\", r"\backslash "), ("&", r"\&"), ("%", r"\%"), ("#", r"\#")):
        s = s.replace(a, b)
    s = _scripts(s)                            # ★прогони індексів — ДО посимвольних
    for ch, m in MATH.items():
        s = s.replace(ch, _m(m))
    for w, m in MATH_WORDS:
        # ★заміна ЛЯМБДОЮ, не рядком: re.sub тлумачить «\mathfrak» у рядку-заміні
        #   як escape-послідовність і падає на bad escape. Пастка тиха: на мапі
        #   без бекслешів вона б не спрацювала й чекала б до першої формули.
        s = re.sub(rf"(?<![A-Za-z\\]){w}(?![A-Za-z])", lambda _mo, _x=m: _m(_x), s)
    return s


def _inline(s):
    """Рядок md → рядок tex. Порядок несучий: спершу код, тоді жирний/курсив."""
    out, i = [], 0
    for m in re.finditer(r"`([^`]+)`", s):
        out.append(_txt(s[i:m.start()]))
        out.append("$" + _math(m.group(1)) + "$")
        i = m.end()
    out.append(_txt(s[i:]))
    r = "".join(out)
    r = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", r)
    r = re.sub(r"(?<![\*\w])\*([^*]+?)\*(?!\*)", r"\\emph{\1}", r)
    return r


# ★Прогін кирилиці → \foreignlanguage{ukrainian}{…}.
#   ЧОМУ ТАК, А НЕ «зробити T2A типовим»: перша збірка повного дала 176 помилок
#   «\cyrn unavailable in encoding T1» — у [T2A,T1] типовим стає ОСТАННЄ, тобто
#   T1, а кириличні гліфи живуть у T2A. Можна було б перевернути порядок, але
#   тоді в T2A поїхав би ВЕСЬ англійський текст (лігатури, переноси). Чесний
#   лік — позначити цитати тим, чим вони є: українським текстом. babel сам
#   перемкне кодування і дасть правильні переноси.
CYR_RUN = re.compile(r"[Ѐ-ӿ]+(?:[\s,.\-!?;:()'`]+[Ѐ-ӿ]+)*")


def _txt(s):
    s = _esc(s)
    for ch, t in TEXT.items():
        s = s.replace(ch, t)
    s = CYR_RUN.sub(lambda mo: r"\foreignlanguage{ukrainian}{" + mo.group(0) + "}", s)
    s = re.sub(r"[⁰-₟¹²³ᵀ]+",
               lambda mo: "$" + _scripts(mo.group(0)) + "$", s)
    for ch, m in MATH.items():                 # математика поза backtick'ами
        s = s.replace(ch, f"${_m(m)}$")
    return s


def gate_glyph(tex, src_name):
    """★ГЕЙТ-ГЛІФ: не-ASCII у ВИХОДІ = символ, якого мапа не знає."""
    bad, cyr = {}, 0
    for i, ln in enumerate(tex.split("\n"), 1):
        for ch in ln:
            if ord(ch) <= 127:
                continue
            if CYRILLIC.match(ch):
                cyr += 1
                continue
            bad.setdefault(ch, i)
    if cyr:
        print(f"  ★ВИНЯТОК {src_name}: кирилиці пропущено {cyr} символів — цитати "
              f"автора ДОСЛІВНО (провенанс). Несе T2A+babel у преамбулі, не викидається.")
    # ★★ГЕЙТ-РОЗМІТКА: залишковий markdown у .tex = розмітка, що НЕ спрацювала.
    #   Чому окремим гейтом: «**» — ASCII, тому гейт-гліф його не бачить;
    #   для LaTeX це валідний текст, тому компілятор мовчить; лог чистий.
    #   Дефект видно ЛИШЕ у PDF — і саме так його й спіймано (витягом тексту).
    #   ⟹ те, що ловиться лише очима, мусить ловитись машиною.
    #   ★ХИБНЯК, знятий на прогоні: перший захід ловив ще й backtick — а «`» у
    #   .tex ЛЕГАЛЬНИЙ: це LaTeX-ські відкривальні лапки, у які я САМА мапаю «“».
    #   Гейт ловив власний легальний вихід (третій такий випадок за день:
    #   fence на «Суд=Омега», export на самому собі, тепер цей). Лишається «**» —
    #   у .tex воно не значить нічого й може приїхати лише як незроблена робота.
    md_left = re.findall(r"\*\*", tex)
    if md_left:
        print(f"  [FAIL] ГЕЙТ-РОЗМІТКА {src_name} — залишкового markdown: "
              f"{len(md_left)} ⟹ у PDF поїхали б літеральні зірочки")
        for mo in list(re.finditer(r".{0,40}\*\*.{0,40}", tex))[:3]:
            print(f"        · «{mo.group(0).strip()}»")
        return len(md_left)

    if bad:
        print(f"  [FAIL] ГЕЙТ-ГЛІФ  {src_name} — не-ASCII у виході: {len(bad)}")
        for ch, i in sorted(bad.items()):
            print(f"        · «{ch}» (U+{ord(ch):04X}) рядок {i} — ДОДАЙ У МАПУ або "
                  f"переформулюй; тихо зникне у PDF")
        return len(bad)
    print(f"  [PASS] ГЕЙТ-ГЛІФ  {src_name} — усі символи мапа знає")
    return 0


PREAMBLE = r"""%% GENERATED from %(src)s -- DO NOT EDIT BY HAND.
%% Edit the source (.md), then run: python src/tools/md_to_tex.py
%% Single column by design: two columns are not convenient for fast reading.
%% NOT COMPILED by the generator: the generating system has no TeX toolchain.
%%   Until a PDF is built somewhere, this file is a claim, not a result.
\documentclass[11pt,a4paper,onecolumn]{article}
%% T2A is required, not decoration: the author's quotations are kept verbatim in
%% Ukrainian (the provenance is his own words), with an English gloss alongside.
\usepackage[T2A,T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage[main=english,ukrainian]{babel}
\usepackage{lmodern}
\usepackage{amsmath,amssymb,amsthm}
\usepackage[margin=2.4cm]{geometry}
\usepackage{microtype}
\usepackage[hidelinks]{hyperref}
\setlength{\parskip}{0.55em}
\setlength{\parindent}{0pt}
%% No \newcommand{\qed}: amsthm defines it, and redefining it is an error.
\title{%(title)s}
\author{%(author)s}
\date{%(date)s}
\begin{document}
\maketitle
"""


def convert(md_path, title, author, date):
    raw = md_path.read_text(encoding="utf-8")
    raw = re.sub(r"<!--.*?-->", "", raw, flags=re.S)      # GUARD-блоки не їдуть у PDF
    body, i, lines = [], 0, raw.split("\n")
    in_list = False
    while i < len(lines):
        ln = lines[i].rstrip()
        if not ln.strip():
            if in_list:
                body.append(r"\end{itemize}")
                in_list = False
            body.append("")
            i += 1
            continue
        if ln.startswith("# ") or ln.strip() == "---":
            i += 1
            continue                                       # тайтл/лінійки — у преамбулі
        if ln.startswith("#### "):
            if in_list:
                body.append(r"\end{itemize}"); in_list = False
            body.append(r"\subsection*{" + _inline(ln[5:]) + "}")
        elif ln.startswith("### "):
            body.append(r"\subsection*{" + _inline(ln[4:]) + "}")
        elif ln.startswith("## "):
            if in_list:
                body.append(r"\end{itemize}"); in_list = False
            body.append(r"\section*{" + _inline(ln[3:]) + "}")
        elif ln.startswith("> "):
            blk = []
            while i < len(lines) and lines[i].startswith(">"):
                blk.append(lines[i].lstrip(">").strip())
                i += 1
            body.append(r"\begin{quote}" + _inline(" ".join(blk)) + r"\end{quote}")
            continue
        elif ln.startswith("- "):
            if not in_list:
                body.append(r"\begin{itemize}"); in_list = True
            item = [ln[2:]]
            i += 1
            while i < len(lines) and lines[i].startswith("  ") and lines[i].strip():
                item.append(lines[i].strip()); i += 1
            body.append(r"\item " + _inline(" ".join(item)))
            continue
        else:
            # ★★АБЗАЦ ЗБИРАЄТЬСЯ ЦІЛИМ, і це не косметика.
            #   Баг, знайдений ВИТЯГОМ ТЕКСТУ З PDF (не логом, не гейтом):
            #   markdown переносить рядки, тому «**consistency, not\n
            #   independent confirmation**» відкриває жирний в одному рядку й
            #   закриває в наступному. Порядкова обробка не матчила його НІКОЛИ
            #   ⟹ у PDF їхали ЛІТЕРАЛЬНІ зірочки. Компілятор мовчав (це валідний
            #   текст), гейт-гліф мовчав (зірочка — ASCII), лог був чистий.
            #   ⟹ та сама ХИБНА ОДИНИЦЯ ВИМІРУ, що й у гейті README: рядок
            #   замість абзацу. Розмітка — властивість АБЗАЦУ.
            para = [ln]
            i += 1
            while i < len(lines):
                nxt = lines[i].rstrip()
                if (not nxt.strip() or nxt.startswith(("#", ">", "- ", "```"))
                        or nxt.strip() == "---"):
                    break
                para.append(nxt)
                i += 1
            body.append(_inline(" ".join(para)))
            continue
        # ★сюди доходять ЛИШЕ гілки заголовків — решта робить continue, просунувши i
        #   сама. Правлячи абзаци, я знесла цей рядок разом зі старим блоком, і «##»
        #   закрутився на місці НАЗАВЖДИ. Спіймано не оком, а тим, що процес не
        #   повернувся: зависання — теж вимір, і краще за тихий неправильний вихід.
        i += 1
    if in_list:
        body.append(r"\end{itemize}")
    # ★title/author НЕ через _inline: там уже готова LaTeX-склейка (\\[0.3em] тощо),
    #   а _inline екранує бекслеші ⟹ перший захід дав у тайтлі
    #   «\textbackslash{}\textbackslash{}[0.3em]» замість переносу. Гейт-гліф цього
    #   НЕ ловить — він міряє СИМВОЛИ, не СТРУКТУРУ; спіймало око на виході.
    head = PREAMBLE % {"src": md_path.name, "title": title,
                       "author": author, "date": date}
    return head + "\n".join(body) + "\n\\end{document}\n"


def main():
    import yaml
    d = yaml.safe_load((EXP / "citation.yaml").read_text(encoding="utf-8"))
    # текст екрануємо, LaTeX-склейку лишаємо сирою
    title = (_txt(d["title"]) + r"\\[0.35em]" + r"\large\normalfont "
             + _txt(f"({d['series']})"))
    author = d["creators"][0]["name"].replace(", ", " ").split()
    author = f"{author[1]} {author[0]}" if len(author) > 1 else author[0]
    date = r"\today" if d.get("date_released") in (None, "TBD") else d["date_released"]

    check = "--check" in sys.argv
    fails = 0
    print("=" * 72)
    print("md_to_tex — одноколонковий .tex з одного джерела")
    print("=" * 72)
    for name in ("compact_en.md", "theory_en.md"):
        src = EXP.parent / name
        if not src.exists():
            print(f"  ⚠ НЕМА {name}")
            continue
        tex = convert(src, title, author, date)
        fails += gate_glyph(tex, name)
        out = EXP / (src.stem + ".tex")
        if fails:
            continue
        if check:
            if not out.exists() or out.read_text(encoding="utf-8") != tex:
                print(f"  [FAIL] {out.name}: розійшовся з {name} (правили руками?)")
                fails += 1
            else:
                print(f"  [PASS] {out.name} ⟷ {name} — збігається")
        else:
            out.write_text(tex, encoding="utf-8")
            print(f"  ✓ {out.name} ({len(tex.splitlines())} рядків)")
    print()
    if fails:
        print(f"✗ НЕ ВИПУЩЕНО: {fails}")
        return 1
    print("⚠★.tex згенеровано, але НЕ СКОМПІЛЬОВАНО: у цій системі немає TeX.")
    print("   Доки PDF не зібрано хоч десь — це КЛЕЙМ, не результат.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
