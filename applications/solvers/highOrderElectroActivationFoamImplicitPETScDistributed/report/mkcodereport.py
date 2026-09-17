#!/usr/bin/env python3
"""Genera `codigo_highOrderManufacturedFDAImplicitPETScDistributed.tex`:
un recorrido del solver, en el orden del archivo, con el codigo a la vista.

Por que el codigo se INCLUYE y no se copia
------------------------------------------
Los listados salen de `\\lstinputlisting` apuntando al `.C` real, con rangos de
lineas. Nada de codigo vive dentro del .tex.

Es la misma decision que en la suite de verificacion: un reporte que copia el
codigo empieza correcto y se vuelve mentira en el primer commit, y peor todavia,
se vuelve mentira EN SILENCIO -- nadie relee un anexo de 60 paginas para
comprobar que sigue coincidiendo. Incluyendo, el PDF de hoy muestra el codigo de
hoy.

Y por eso los rangos NO estan escritos a mano: este script BUSCA cada funcion por
nombre y calcula su rango contando llaves. Si alguien agrega cincuenta lineas
arriba, se regenera y los rangos se recalculan solos. Un rango a mano seria
exactamente el mismo problema que copiar el codigo, disfrazado.

Uso
---
    python3 mkcodereport.py && pdflatex codigo_*.tex
"""

from __future__ import annotations

import datetime as _dt
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOLVER = HERE.parent
NAME = SOLVER.name
SOURCE = NAME + ".C"
FIELDS = "createFields.H"
OUT = HERE / f"codigo_{NAME}.tex"

MAX_LISTING = 60          # lineas por listado antes de recortar

# `listings` es orientado a BYTES: un caracter UTF-8 multibyte dentro de un
# listado rompe la compilacion con "Invalid UTF-8 byte sequence", y el mensaje no
# dice que archivo ni que caracter. Se mapean uno por uno con `literate`.
#
# La tabla se VERIFICA contra el fuente al generar: si aparece un caracter no
# ASCII que no este aca, el script aborta nombrandolo. Filtrarlo en silencio
# dejaria un listado con un caracter comido, que es peor que no compilar --
# nadie revisa un anexo de sesenta paginas buscando un simbolo faltante.
LITERATE = {
    "\u00d7": r"{$\times$}",
    "\u2248": r"{$\approx$}",
    "\u2014": r"{---}",
    "\u2192": r"{$\rightarrow$}",
    "\u2264": r"{$\leq$}",
    "\u2265": r"{$\geq$}",
    "\u00b0": r"{$^\circ$}",
}


def literate_clause(text: str) -> str:
    """La clausula `literate=` de listings, y un chequeo de cobertura."""
    unknown = sorted({c for c in text if ord(c) > 127} - set(LITERATE))
    if unknown:
        raise SystemExit(
            "El fuente trae caracteres no ASCII sin mapear en LITERATE: "
            + ", ".join(f"{c!r} (U+{ord(c):04X})" for c in unknown)
            + ". Agregalos a la tabla; listings no los puede imprimir solo.")
    used = [c for c in LITERATE if c in text]
    if not used:
        return ""
    return "  literate=" + " ".join(f"{{{c}}}{LITERATE[c]}1" for c in used) + ",\n"


# ---------------------------------------------------------------------------
# localizar funciones en el fuente
# ---------------------------------------------------------------------------

def find_functions(path: Path) -> dict:
    """{nombre: (primera_linea, ultima_linea)} contando llaves.

    Acepta las dos formas que el archivo usa: `tipo nombre(args)` en una linea y
    `tipo nombre` con el parentesis en la siguiente, que es la mas comun aca.
    """
    lines = path.read_text().split("\n")
    out = {}
    for i in range(len(lines) - 1):
        line = lines[i]
        if line.strip().startswith("//"):
            continue
        m = re.match(r"^(\s*)(?:static\s+|inline\s+)?"
                     r"([A-Za-z_][\w:<>,&\* ]*?)\s+([A-Za-z_]\w*)\s*$", line)
        if not (m and re.match(r"^\s*\(\s*$", lines[i + 1])):
            m = re.match(r"^(\s*)(?:static\s+|inline\s+)?"
                         r"([A-Za-z_][\w:<>,&\* ]*?)\s+([A-Za-z_]\w*)\s*\([^;]*\)\s*$",
                         line)
            if not m:
                continue
        name = m.group(3)
        # el cuerpo arranca en la primera `{` a la altura de la definicion
        depth, start_body = 0, None
        for j in range(i, min(i + 80, len(lines))):
            if "{" in lines[j]:
                start_body = j
                break
            if ";" in lines[j] and "(" not in lines[j]:
                break
        if start_body is None:
            continue
        for j in range(start_body, len(lines)):
            depth += lines[j].count("{") - lines[j].count("}")
            if depth == 0 and j > start_body - 1:
                out.setdefault(name, (i + 1, j + 1))
                break
    return out


def span(funcs, name, limit=MAX_LISTING):
    """(primera, ultima, recortada?) del rango a listar."""
    if name not in funcs:
        raise KeyError(f"no encontre `{name}` en el fuente")
    a, b = funcs[name]
    if b - a + 1 > limit:
        return a, a + limit - 1, True
    return a, b, False


# ---------------------------------------------------------------------------
# el recorrido
# ---------------------------------------------------------------------------
# (titulo, prosa, [(nombre_funcion | (primera,ultima), pie)])

SECTIONS = [
 ("Qué resuelve este solver", r"""
Resuelve la propagación del potencial de acción cardíaco: la ecuación de
monodominio acoplada al modelo iónico \textbf{ten Tusscher--Noble--Noble--Panfilov
(TNNP)}, con el objetivo de reproducir el benchmark de Niederer et al. (2012).

A diferencia del solver MMS, acá \textbf{no hay solución exacta}. La verificación
es por convergencia, por invariancia frente a la descomposición, y contra
referencia publicada.

En cada paso de tiempo hay dos trabajos muy distintos:
%
\begin{enumerate}
  \item integrar, en cada celda o punto de Gauss, un sistema \textbf{rígido} de
        19 EDO --- el modelo TNNP --- para obtener la corriente iónica $I_{ion}$;
  \item resolver el sistema lineal implícito de la difusión.
\end{enumerate}
%
El primero domina el tiempo: la resolución lineal completa es apenas el 17--44\,\%
del lazo. Eso explica varias decisiones del archivo, entre ellas que la única
región OpenMP esté sobre las EDO y no sobre el ensamblado.

\textbf{Dos ejes de paralelismo}, en lugares distintos del archivo: \textbf{MPI}
reparte la malla entre ranks y toca todas las fases; \textbf{OpenMP} reparte el
lazo de puntos de integración dentro de un rank y sólo toca las EDO.
""", []),

 ("El modelo iónico TNNP", r"""
El modelo vive en una clase dentro del archivo, con sus 19 estados y sus
variables algebraicas.

\texttt{computeRates} evalúa las derivadas; \texttt{computeVariables} las
cantidades algebraicas; \texttt{calculateCurrent} devuelve $I_{ion}$.

\texttt{protectState} es la protección numérica: acota compuertas y
concentraciones a rangos físicos. Lleva un contador de correcciones, y ése es el
único estado compartido de todo el camino de las EDO --- por eso su incremento
va bajo \texttt{omp critical} y no \texttt{atomic}: el test de cota y el
incremento tienen que ser un solo paso indivisible.

El test \textbf{E01} verifica este modelo aislado, en 0-D, contra los valores
publicados de ten Tusscher 2004 --- APD90, pico, reposo, dV/dt máximo. Es el
test más barato de la suite y sostiene todo lo demás: si la corriente iónica
estuviera mal, ningún test de convergencia lo detectaría, porque todos
convergerían limpiamente a la respuesta equivocada.
""", [("protectState", "La protección numérica y su contador compartido."),
      ("calculateCurrent", "La corriente iónica. Recortada.")]),

 ("Los integradores de EDO, y el único paralelismo por hilos", r"""
El solver trae sus propios steppers en vez de usar el marco \texttt{ODESolver}
de OpenFOAM: \texttt{eulerStep}, \texttt{rk4Step} y \texttt{rkf45Trial}
adaptativo. Un nombre desconocido es \texttt{FatalError}, no un fallback
silencioso.

Euler y RK4 dan \textbf{un} paso del tamaño del paso exterior completo, así que
su exactitud la controla \texttt{dt}; RKF45 subdivide por su cuenta y la
controla su tolerancia. El test \textbf{E02} hace converger los tres al mismo
potencial de acción, que es una afirmación mucho más fuerte que apretar uno
solo contra sí mismo.

El paso adaptativo \textbf{persiste entre llamadas} en \texttt{stepMs\_}, de modo
que cada llamada hereda el $h$ de la anterior. Eso hacía que la clave
\texttt{initialODEStep} fuera inalcanzable --- el constructor siembra
\texttt{stepMs\_} positivo, así que la rama que la leía sólo podía dispararse con
\texttt{deltaT == 0}. El test \textbf{E03} lo midió (cero efecto sobre cuatro
órdenes de magnitud) y la clave se removió.

\texttt{solveODE} contiene la \textbf{única} región OpenMP del archivo. Los
hilos reparten los puntos de integración, y funciona porque cada punto resuelve
su propia EDO sin hablar con ningún otro --- por eso el resultado es
bit-idéntico con cualquier cantidad de hilos, cosa que \textbf{E20} mide dando
exactamente cero. \textbf{El ensamblado no está paralelizado con hilos.}
""", [("rkf45Trial", "El paso adaptativo. Recortado."),
      ("solveODE", "El despacho de stepper y la región OpenMP. Recortado.")]),

 ("El residuo de estados: una colectiva que colgaba", r"""
\texttt{maxStateRelativeL2Difference} da el cambio relativo máximo de los
estados entre dos iteraciones no lineales, y es parte del criterio de parada.

Es una función \textbf{colectiva}: reduce dos veces por estado. Y ahí había un
cuelgue. La versión anterior devolvía temprano cuando el arreglo de estados de un
rank venía vacío --- exactamente lo que le pasa a un rank cuya partición no tiene
celdas --- así que ese rank salía sin participar de las reducciones mientras los
demás quedaban bloqueados esperando. \textbf{No es un número mal calculado: es
que la corrida no termina.}

Está corregido de dos formas, porque cualquiera de las dos sola dejaba el cuelgue
vivo: el rank vacío ahora cae hasta los lazos de acumulación (que iteran cero
veces) y llega igual a los dos \texttt{reduce}, y el conteo de vueltas del lazo
se reduce con \texttt{minOp}, porque salía de \texttt{current[0].size()}, que es
rank-local y ni siquiera se puede leer en un rank vacío.

El test \textbf{E19} lo mide sobre una malla de dos celdas: antes moría en el
timeout de 300 s a np=4 y np=8, ahora corre en 1 s, igual que el control serial.
""", [("maxStateRelativeL2Difference", "La colectiva, con las dos correcciones.")]),

 ("Bloques de construcción del ensamblado", r"""
Antes de las matrices de rigidez vienen las piezas compartidas.

\texttt{addTripletIfNeeded} es el punto por el que pasa toda entrada de la
matriz, y descarta las menores que \texttt{SMALL}.

\texttt{addCellGradientDotCoeffs} agrega la fila de gradiente de una celda
proyectada sobre un vector, usando los stencils MLS. Es donde vive el problema
abierto más grande: sobre celdas pegadas a un corte de procesador ese stencil
pierde un miembro, y como el ajuste se rehace sobre el conjunto reducido,
\textbf{refita toda la fila}. Vive en solids4foam, no acá.

\texttt{exchangeCoupledFaceRows} intercambia con el rank vecino las filas que el
término de estabilización necesita del otro lado del corte. Es \textbf{colectiva}
y está deliberadamente \emph{fuera} del lazo de chunks: adentro se trabaría en
cuanto dos ranks tuvieran distinta cantidad de chunks.
""", [("addTripletIfNeeded", "Por acá pasa toda entrada."),
      ("exchangeCoupledFaceRows", "El intercambio a través del corte. Recortado.")]),

 ("Las dos matrices de rigidez, y el defecto de la conductividad", r"""
Igual que en el MMS hay dos ensamblados y la tríada elige cuál corre: el
ortogonal de dos puntos y el de alto orden sobre el stencil LRE.

Acá estuvo el defecto más grande que encontró la verificación. El camino de alto
orden usaba \texttt{conductivity[own]} \textbf{en las dos filas} de una cara
interna, en vez del promedio $\tfrac12(D_{own}+D_{nei})$ que ya usaban el camino
de bajo orden y el término de estabilización.

Con $D$ uniforme los dos tensores coinciden y la diferencia se cancela
idénticamente, así que fue invisible durante toda la vida del código: todas las
corridas usaban el tensor constante del benchmark. Con $D$ \textbf{no uniforme}
--- o sea, en cuanto entren fibras --- hacía que $K$ cambiara un \textbf{21\,\%}
al descomponer la malla, porque al cortar, cada rank pasa a ser dueño de su lado
de la cara y ninguno coincide con lo que hizo la corrida serial.

El arreglo son tres cambios, no uno: el promedio en la cara interna, el promedio
\textbf{también en la cara de procesador} (con el tensor remoto), y sacar el
intercambio de conductividad de la condición de estabilización, porque ahora hace
falta con cualquier alpha. Corregir sólo el primero arreglaba el operador serial
y dejaba el paralelo roto. \textbf{M15} lo mide: de 2.108e-01 a 4.736e-15.
""", [("assembleHighOrderStiffnessMatrix", "El ensamblado de alto orden. Recortado.")]),

 ("Modo híbrido: EDO donde el frente lo pide", r"""
Resolver la EDO en cada punto de Gauss es lo más directo y lo más caro;
resolverla en el centro de celda y reconstruir a los puntos es mucho más barato
pero falla donde el frente es abrupto, porque los estados del TNNP son rígidos.

\texttt{flagFrontCells} clasifica las celdas del frente por el rango de $V_m$
entre sus puntos de Gauss, y \texttt{dilateFrontCells} agranda esa banda unas
capas. El modo híbrido resuelve la EDO en los puntos de Gauss ahí y en el centro
de celda en el resto.

El detalle de implementación que importa: la clasificación se calcula
\textbf{una vez por paso} y se congela para todo el solve no lineal. Recalcularla
por iteración hace que el interruptor frente/suave oscile entre iteraciones y
convierte el punto fijo de Picard en un ciclo límite que nunca converge.
""", [("flagFrontCells", "La clasificación del frente. Recortada."),
      ("dilateFrontCells", "La dilatación de la banda. Recortada.")]),

 ("Tiempos de activación y muestreo del benchmark", r"""
Lo que el benchmark reporta son \textbf{tiempos de activación}, no campos.

\texttt{updateActivationTimes} los interpola \emph{dentro} del paso ---
$t_n + w\,\Delta t$ con $w$ el cruce lineal del umbral --- así que no quedan
cuantizados a múltiplos de $\Delta t$, lo cual importa para cualquier estudio de
refinamiento.

\texttt{writeActivationSamples} escribe el valor en el punto P8 y un perfil a lo
largo de la diagonal. El muestreo es en puntos \textbf{físicos fijos},
independientes de la malla, y eso es lo que permite restar dos perfiles índice a
índice sin interpolar entre discretizaciones --- la base de E11 y de los estudios
de convergencia.

\texttt{nearestCellToPoint} merece atención: devuelve el índice en \textbf{exactamente
un} rank y $-1$ en todos los demás, resuelto con una reducción para que la
elección sea determinista. El consumidor está guardado contra ese $-1$; el
\emph{banner} que imprime el punto de parada no lo estaba, y evaluaba
\texttt{mesh.C()[-1]} en todos los ranks perdedores. En build optimizado eso es
lectura fuera de rango, o sea comportamiento indefinido, y se manifestaba como
\textbf{SIGFPE} de forma reproducible pero aparentemente arbitraria: a dx=0.125
mm fallaba con 4 ranks y no con 1, 2 ni 8.
""", [("updateActivationTimes", "La interpolación dentro del paso. Recortada."),
      ("nearestCellToPoint", "Exactamente un rank gana; el resto devuelve -1.")]),

 ("main(): montaje y lazo temporal", r"""
\texttt{main()} lee malla y diccionarios, construye los operadores LRE, ensambla
$M$ y $K$ \textbf{una sola vez} --- la no linealidad está en el término iónico,
no en la difusión --- y entra al lazo temporal.

Cada paso resuelve un sistema no lineal con Picard, JFNK o diagonalIion, y
\textbf{E10} verifica que los tres lleguen al mismo lugar. Un nombre de método
inválido es \texttt{FatalError} y no un fallback silencioso a Picard, cosa que
\textbf{E18} comprueba.

Al agotar las iteraciones sin converger, la rama está \textbf{una sola vez},
después del despacho de método, así que los tres comparten la semántica: por
defecto se descarta el paso y se vuelve a los valores del inicio, y con
\texttt{nonlinearAcceptUnconverged} se acepta el iterado. \textbf{E24} lo
verifica en las dos configuraciones. El solver MMS tiene el corte adentro del
despacho y por eso allá Picard acepta mientras los otros dos descartan --- dos
semánticas en la misma matriz de métodos, y sigue abierto.
""", [("computeStableDeltaT", "El paso estable de referencia."),
      ("applyStimulus", "El estímulo a nivel tejido. Recortado.")]),
]


# ---------------------------------------------------------------------------
# emision
# ---------------------------------------------------------------------------

PREAMBLE = r"""\documentclass[11pt,a4paper]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[spanish,es-noquoting,es-noshorthands]{babel}
\usepackage{amsmath}
\usepackage[margin=2.2cm]{geometry}
\usepackage{listings}
\usepackage{xcolor}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage[hidelinks]{hyperref}
\setlength{\parskip}{0.55em}
\setlength{\parindent}{0pt}

\definecolor{cmt}{rgb}{0.35,0.50,0.35}
\definecolor{kw}{rgb}{0.10,0.25,0.60}
\definecolor{str}{rgb}{0.60,0.20,0.20}
\definecolor{bg}{rgb}{0.975,0.975,0.97}

\lstset{
  language=C++,
  basicstyle=\ttfamily\scriptsize,
  commentstyle=\color{cmt}\itshape,
  keywordstyle=\color{kw}\bfseries,
  stringstyle=\color{str},
  backgroundcolor=\color{bg},
  numbers=left, numberstyle=\tiny\color{gray}, numbersep=7pt,
  breaklines=true, breakatwhitespace=false,
  showstringspaces=false, tabsize=4,
  frame=leftline, framesep=6pt, rulecolor=\color{gray},
  columns=fullflexible, keepspaces=true,
  captionpos=b,
LITERATE_HERE}
"""


def emit_listing(rel_source, first, last, caption, truncated):
    cap = caption + (r" \textit{(recortado; sigue en el fuente.)}" if truncated else "")
    return "\n".join([
        r"\begin{lstlisting}[firstnumber=%d, firstline=%d, lastline=%d, caption={%s}]"
        % (first, first, last, cap),
        r"\end{lstlisting}",
    ]).replace(r"\begin{lstlisting}",
               r"\lstinputlisting[firstnumber=%d, firstline=%d, lastline=%d, caption={%s}]{%s}%%"
               % (first, first, last, cap, rel_source)).replace(r"\end{lstlisting}", "")


def main():
    src = SOLVER / SOURCE
    text = src.read_text()
    funcs = find_functions(src)
    rel = "../" + SOURCE

    body = []
    for title, prose, items in SECTIONS:
        body.append(r"\section{%s}" % title)
        body.append(prose.strip())
        for ref, caption in items:
            if isinstance(ref, tuple):
                a, b, trunc = ref[0], ref[1], False
            else:
                try:
                    a, b, trunc = span(funcs, ref)
                except KeyError as exc:
                    print(f"  aviso: {exc}", file=sys.stderr)
                    continue
            body.append(emit_listing(rel, a, b, caption, trunc))

    # indice completo de funciones
    index = [r"\section{Índice de funciones}",
             "Todas las funciones del espacio de nombres anónimo, en orden de "
             "aparición. Las secciones anteriores recorren las principales; ésta "
             "es el mapa completo para ubicarse en el archivo.",
             r"\begin{longtable}{rll}", r"\toprule",
             r"línea & función & líneas \\", r"\midrule", r"\endhead"]
    for name, (a, b) in sorted(funcs.items(), key=lambda kv: kv[1][0]):
        index.append(r"%d & \texttt{%s} & %d \\" % (a, name.replace("_", r"\_"), b - a + 1))
    index += [r"\bottomrule", r"\end{longtable}"]

    tex = "\n\n".join([
        PREAMBLE.replace('LITERATE_HERE', literate_clause(text)),
        r"\title{\textbf{Recorrido del código}\\[4pt]\large\texttt{%s}}"
        % NAME.replace("_", r"\_"),
        r"\author{}",
        r"\date{Generado el %s a partir de \texttt{%s} (%d líneas)}"
        % (_dt.date.today().isoformat(), SOURCE.replace("_", r"\_"),
           len(src.read_text().split("\n"))),
        r"\begin{document}", r"\maketitle",
        r"\noindent\textbf{Cómo leer esto.} Los listados NO son copias: salen de "
        r"\texttt{\textbackslash lstinputlisting} apuntando al fuente real, con "
        r"rangos que este reporte recalcula al regenerarse. El PDF de hoy muestra "
        r"el código de hoy. Los números de línea del margen son los del archivo.",
        r"\tableofcontents", r"\newpage",
        "\n\n".join(body),
        "\n".join(index),
        r"\end{document}",
    ])
    OUT.write_text(tex)
    print(f"{OUT.name}: {len(SECTIONS)} secciones, {len(funcs)} funciones indexadas")
    return 0


if __name__ == "__main__":
    sys.exit(main())
