#!/usr/bin/env python3
"""Genera `codigo_<solver>.tex`: un recorrido del solver con el codigo a la vista.

Desde 2026-10 buena parte del codigo vive en `src/monodomainCommon/` (headers
compartidos con el solver MMS, incluidos dentro del espacio de nombres anonimo
del .C). El script busca cada funcion primero en el .C y despues en esos
headers, y cada listado apunta al archivo donde la funcion vive.

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
COMMON = SOLVER.parent.parent.parent / "src" / "monodomainCommon"

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

def find_functions(path: Path, start: int = 0) -> dict:
    """{nombre: (primera_linea, ultima_linea)} contando llaves.

    Acepta las dos formas que el archivo usa: `tipo nombre(args)` en una linea y
    `tipo nombre` con el parentesis en la siguiente, que es la mas comun aca.
    """
    lines = path.read_text().split("\n")
    out = {}
    for i in range(start, len(lines) - 1):
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
 ("Qué resuelve este solver, y dónde vive el código", r"""
Resuelve la propagación del potencial de acción cardíaco: la ecuación de
monodominio acoplada al modelo iónico \textbf{ten Tusscher--Noble--Noble--Panfilov
(TNNP)}, con el objetivo de reproducir el benchmark de Niederer et al. (2012).
A diferencia del solver MMS, acá \textbf{no hay solución exacta}: la verificación
es por convergencia de Cauchy, por invariancia frente a la descomposición y
contra referencia publicada.

En cada paso de tiempo hay dos trabajos muy distintos:
%
\begin{enumerate}
  \item integrar, en cada celda o punto de Gauss, un sistema \textbf{rígido} de
        19 EDO --- el modelo TNNP --- para obtener la corriente iónica $I_{ion}$;
  \item resolver el sistema lineal implícito de la difusión.
\end{enumerate}
%
El primero domina el tiempo. \textbf{MPI} reparte la malla entre ranks y toca
todas las fases; \textbf{OpenMP} reparte, dentro de un rank, el lazo de EDO por
punto y las reconstrucciones en los puntos de Gauss (de $V_m$ y de los
estados), además de los productos matriz--vector de Eigen.

\textbf{Dónde vive el código (2026-10).} Todo lo que el solver comparte con el MMS
--- álgebra distribuida y PETSc, ensamblado de masa y rigidez, reconstrucciones,
GMRES/JFNK, esquemas temporales, criterio de parada --- está en
\texttt{src/monodomainCommon/}, en headers que se incluyen \emph{dentro} del
espacio de nombres anónimo del \texttt{.C}: no hay biblioteca que enlazar, y un
cambio ahí cambia los dos solvers. En el \texttt{.C} queda lo propio: el modelo
TNNP, el estímulo, el modo híbrido, los tiempos de activación y el lazo temporal.
Cada listado de este recorrido indica de qué archivo sale.
""", []),

 ("El modelo iónico TNNP", r"""
El modelo vive en una clase dentro del archivo, con sus 19 estados y sus
variables algebraicas. \texttt{computeRates} evalúa las derivadas;
\texttt{computeVariables} las cantidades algebraicas; \texttt{calculateCurrent}
devuelve $I_{ion}$.

\texttt{protectState} es la protección numérica: acota compuertas a $[0,1]$,
pone pisos a las concentraciones y a $V_m$. Lleva un contador de correcciones,
el único estado compartido de todo el camino de las EDO --- por eso su
incremento va bajo \texttt{omp critical}. \textbf{Está apagada por defecto}, en
el solver y desde 2026-10-16 también en el harness: encendida, su contador llega
al tope de $10^6$ y los BDF pierden el orden (E02: con estados CC, BDF2/3/4 dan
1.74/2.22/3.54 con ella y 1.89/2.96/4.60 sin ella), y ninguna corrida la
necesitó.

El test \textbf{E26} verifica el modelo aislado, en 0-D, contra los valores
publicados de ten Tusscher 2004 (el pico, contra el propio modelo convergido).
""", [("protectState", "La protección numérica y su contador compartido."),
      ("calculateCurrent", "La corriente iónica. Recortada.")]),

 ("Los integradores de EDO y el driver de Vm", r"""
El solver trae sus propios steppers: \texttt{eulerStep}, \texttt{rk4Step} y
\texttt{rkf45Trial} adaptativo. Un nombre desconocido es \texttt{FatalError}.
Cada integración RKF45 arranca desde \texttt{stateODEInitialStep} (el viejo
\texttt{stepMs\_} que heredaba el paso entre llamadas, y hacía inalcanzable la
clave, ya no existe). \textbf{E27} hace converger los tres integradores al mismo
potencial de acción.

Dentro del paso, la EDO no ve un $V_m$ congelado: \texttt{computeRates} lo toma
del \emph{driver} (\texttt{drivenVmMilliVolt}), un polinomio en el tiempo que
pasa por $V^n$ y el candidato $V^{n+1}$ --- y, con \texttt{stateODEDriver
history}, también por los niveles anteriores de la historia BDF o por los valores
de etapa ESDIRK. Con $V_m$ congelado el TNNP daba orden 1 incluso con
Crank--Nicolson; con la recta, 2; la historia es la que deja a BDF3 y BDF4 dar su
orden.

\texttt{solveODE} reparte los puntos con OpenMP; cada punto es una EDO
independiente, así que el resultado es bit a bit igual con cualquier cantidad de
hilos (\textbf{E33}).
""", [("rkf45Trial", "El paso adaptativo. Recortado."),
      ("drivenVmMilliVolt", "El driver de Vm dentro del paso."),
      ("interpolateVmDriver", "La interpolación de Lagrange del driver (compartida)."),
      ("solveODE", "El despacho de stepper y la región OpenMP. Recortado.")]),

 ("Esquemas temporales: theta, BDF1-4 y el arranque ESDIRK3", r"""
\texttt{timeSchemeCoeffs} da los coeficientes de backwardEuler, crankNicolson
y BDF1--BDF4. Un BDF-$k$ necesita $k-1$ niveles previos: los da un arranque
\textbf{ESDIRK3(2)4L[2]SA} (\texttt{esdirk3Tableau}), de tercer orden respecto del
sistema semidiscreto, que es lo que conserva el orden nominal.

La historia se \textbf{reinicia} (y los $k-1$ pasos siguientes vuelven a ser
ESDIRK3) cuando cambia el estímulo --- un salto en el forzado que un método
multipaso no puede cruzar --- en un paso revertido, y en un paso de longitud no
nominal. El estímulo se promedia sobre el paso (\texttt{applyStimulus}), así que
la carga inyectada es exacta para cualquier $\Delta t$; antes duraba un paso de
más.

\textbf{E02} mide el orden de los seis esquemas en cuatro caminos de código (bajo
orden, CC, híbrido y GP). En el par más fino, CN da 2.00 y BDF4 entre 4.6 y 5.6.
""", [("timeSchemeCoeffs", "Coeficientes de los esquemas (compartido)."),
      ("esdirk3Tableau", "La tabla del arranque ESDIRK3 (compartido)."),
      ("applyStimulus", "El estímulo promediado sobre el paso. Recortado.")]),

 ("Reconstrucciones en los puntos de Gauss", r"""
El término fuente es de alto orden porque $I_{ion}$ se evalúa en los puntos de
Gauss y se integra con cuadratura: hace falta $V_m$ (y, en el modo CC, los
estados) en esos puntos. \texttt{reconstructCellFieldsAtIntegrationPoints} lo
hace para \textbf{K campos a la vez}, recorriendo el estencil MLS una sola vez por
celda: cada vecino se ubica una vez (local o remoto) y los coeficientes de
gradiente, hessiano y tercera derivada se cargan una vez para todos los campos.
La primera versión compartida hacía un campo por llamada y dejó la
reconstrucción de estados del MMS un 30\,\% más lenta; para el TNNP, con 19
estados, la ganancia fue mayor. El resultado es bit a bit igual a reconstruir
de a uno.

\texttt{reconstructVmLevelAtIionIntegrationPoints} lleva a los puntos de Gauss
un nivel guardado de $V_m$ (de la historia BDF o de una etapa ESDIRK). Es lo que
usa \texttt{stateODEGaussPointHistory}: con ella, los caminos de puntos de Gauss
(GP y el frente del híbrido) usan el mismo polinomio de la historia que el
camino por celda, en vez de la recta entre $V^n$ y el candidato. Sin ella BDF3 y
BDF4 quedaban en $\sim$2.4 en esos caminos; con ella dan su orden (en el MMS, con
solución exacta, BDF4 3.89/3.91), por un 6--9\,\% de tiempo (\textbf{E36}).
""", [("reconstructCellFieldsAtIntegrationPoints", "El núcleo de K campos (compartido). Recortado."),
      ("reconstructVmLevelAtIionIntegrationPoints", "Un nivel de la historia en los puntos de Gauss (compartido)."),
      ("reconstructStatesAtIionIntegrationPoints", "Los estados del TNNP, todos en una pasada.")]),

 ("El residuo de estados: una colectiva que colgaba", r"""
\texttt{maxStateRelativeL2Difference} da el cambio relativo máximo de los
estados entre dos iteraciones no lineales. Es \textbf{colectiva}: reduce dos
veces por estado. La versión anterior devolvía temprano cuando el arreglo de
estados de un rank venía vacío --- lo que le pasa a un rank sin celdas --- y ese
rank salía sin participar de las reducciones: la corrida no terminaba. Ahora el
rank vacío llega igual a los dos \texttt{reduce}, y el conteo de vueltas se
reduce con \texttt{minOp}. Lo mide \textbf{E32}.
""", [("maxStateRelativeL2Difference", "La colectiva, con las dos correcciones.")]),

 ("El criterio de parada no lineal, común a todos los métodos", r"""
Picard, diagonalIion, JFNK y las etapas ESDIRK paran con el mismo test,
\texttt{nonlinearIterationExit}: residuo acoplado de la ecuación del paso más los
incrementos de $V_m$, $I_{ion}$ y, si se pide, de los estados. El MMS paraba
Picard sólo por incrementos; un test para todos hace que una comparación de
métodos (\textbf{E08}) compare métodos y no reglas de parada.

Tiene una salida por \textbf{piso}: si los incrementos ya pasan y el residuo
acoplado dejó de bajar, la iteración se corta y se informa. Se agregó cuando, con
el criterio común, Picard del MMS repetía 2000 veces la misma solución lineal
con el residuo clavado en $3.66\cdot10^{-12}$ contra $10^{-12}$. La causa era
otra --- el KSP ensamblado paraba sobre el residuo \emph{precondicionado}, y ahora
precondiciona por derecha --- pero el piso queda como red de seguridad.
""", [("nonlinearIterationExit", "El test de parada y la salida por piso (compartido)."),
      ("nonlinearIterationStops", "El conteo y el aviso de una salida por piso (compartido).")]),

 ("Bloques de construcción del ensamblado", r"""
Las piezas compartidas del ensamblado están en \texttt{monodomainAssembly.H} y
\texttt{monodomainBase.H}. \texttt{addTripletIfNeeded} es el punto por el que
pasa toda entrada de la matriz y descarta las menores que \texttt{SMALL}.

\texttt{exchangeCoupledFaceRows} intercambia con el rank vecino las filas que el
término de estabilización necesita del otro lado del corte. Es \textbf{colectiva}
y está deliberadamente \emph{fuera} del lazo de chunks: adentro se trabaría en
cuanto dos ranks tuvieran distinta cantidad de chunks.
""", [("addTripletIfNeeded", "Por acá pasa toda entrada (compartido)."),
      ("exchangeCoupledFaceRows", "El intercambio a través del corte (compartido). Recortado.")]),

 ("Las dos matrices de rigidez, y el defecto de la conductividad", r"""
Hay dos ensamblados y la tríada elige cuál corre: el ortogonal de dos puntos y el
de alto orden sobre el estencil LRE, los dos compartidos con el MMS.

El camino de alto orden usaba \texttt{conductivity[own]} en las dos filas de una
cara interna, en vez del promedio $\tfrac12(D_{own}+D_{nei})$. Con $D$ uniforme
la diferencia se cancela; con $D$ no uniforme hacía que $K$ cambiara un
\textbf{21\,\%} al descomponer la malla. El arreglo fueron tres cambios: el
promedio en la cara interna, el promedio también en la de procesador, y sacar el
intercambio de conductividad de la condición de estabilización. \textbf{M15} lo
mide: de 2.108e-01 a 4.736e-15.
""", [("assembleHighOrderStiffnessMatrix", "El ensamblado de alto orden (compartido). Recortado.")]),

 ("La solución lineal: PETSc con precondicionado por derecha", r"""
\texttt{PetscKspMatrixSolver} guarda la matriz, el KSP y la factorización entre
pasos, y decide en cada solve si reconstruir, refrescar valores o reusar. En
paralelo, \texttt{ilu}/\texttt{lu} se vuelven block-Jacobi con esa factorización
por bloque. Desde 2026-10, GMRES y BiCGStab precondicionan \textbf{por derecha}
(y cg monitorea la norma sin precondicionar): así la tolerancia lineal mide el
mismo residuo verdadero que el criterio no lineal.
""", [(("monodomainPetsc.H", "@// Stop on the TRUE residual", "setPetscOptionsPrefix(ksp_, optionsPrefix);"), "El lado del precondicionador en PetscKspMatrixSolver::reset (compartido).")]),

 ("Modo híbrido: GP en el frente, CC en el resto", r"""
\texttt{flagFrontCells} marca como frente las celdas cuyo $V_m$ varía más de
\texttt{hybridFrontVmRange} entre sus puntos de Gauss (o que están estimuladas),
y \texttt{dilateFrontCells} agranda la banda unas capas. La clasificación se
calcula \textbf{una vez por paso} y se congela para todo el solve no lineal:
recalcularla por iteración convierte el punto fijo en un ciclo límite.

En las celdas del frente la EDO se integra en cada punto de Gauss; en el resto,
una vez en el centro, como en CC puro. \textbf{Desde 2026-10-16 las celdas suaves
guardan su propio estado en el centro} (\texttt{hybridCellStates}) y sólo el
frente vive en los puntos de Gauss. Antes, cada paso rehacían su estado como el
promedio de su reconstrucción en Gauss --- que no es la identidad: agrega
$(h^2/24)\,\Delta s$ --- y eso era una difusión de los estados aplicada una vez por
paso, más fuerte cuanto más chico $\Delta t$: un piso de $2.2\cdot10^{-3}$ en la
convergencia temporal que ningún esquema podía bajar. Con el arreglo el híbrido
da lo mismo que GP con historia (CN orden 2.00, BDF4 5.8e-6 en el par fino) a la
mitad del costo, sin limitador.

Al final de cada corrida con híbrido el log informa cuántas veces cambiaron de
modo las celdas.
""", [("flagFrontCells", "La clasificación del frente. Recortada."),
      ("dilateFrontCells", "La dilatación de la banda. Recortada."),
      ("averageStatesIntegrationPointsToCells", "El promedio de Gauss a celda, que ahora se usa una sola vez por celda.")]),

 ("Tiempos de activación y muestreo del benchmark", r"""
Lo que el benchmark reporta son \textbf{tiempos de activación}.
\texttt{updateActivationTimes} los interpola \emph{dentro} del paso ---
$t_n + w\,\Delta t$ con $w$ el cruce lineal del umbral ---, así que no quedan
cuantizados a múltiplos de $\Delta t$; esa interpolación tiene su propio error
$O(\Delta t^2)$, y por eso E02 mide órdenes sobre el campo $V_m$ y no sobre los
tiempos.

\texttt{writeActivationSamples} escribe el valor en P8 y un perfil sobre la
diagonal, en puntos físicos fijos muestreados por distancia inversa.
\texttt{nearestCellToPoint} devuelve el índice en \textbf{exactamente un} rank y
$-1$ en el resto; el banner que imprimía el punto de parada no estaba guardado
contra ese $-1$ y daba un SIGFPE aparentemente arbitrario. \textbf{E28} verifica
que el frente cruce los cortes de procesador sin déficit.
""", [("updateActivationTimes", "La interpolación dentro del paso. Recortada."),
      ("nearestCellToPoint", "Exactamente un rank gana; el resto devuelve -1.")]),

 ("main(): montaje y lazo temporal", r"""
\texttt{main()} lee malla y diccionarios, construye los operadores LRE, ensambla
$M$ y $K$ una sola vez y entra al lazo temporal. Cada paso resuelve un sistema no
lineal con Picard, JFNK o diagonalIion; \textbf{E08} verifica que los tres lleguen
al mismo lugar, y un nombre inválido es \texttt{FatalError} (\textbf{E31}).

Al agotar las iteraciones sin converger, los tres métodos --- en los dos solvers
--- comparten la semántica: por defecto se descarta el paso y se vuelve a los
valores del inicio (y se reinicia la historia BDF), y con
\texttt{nonlinearAcceptUnconverged} se acepta el iterado. Lo verifica
\textbf{E34}. Los valores por defecto del harness desde 2026-10-16: sin
protección del TNNP, sin limitador de estados, con
\texttt{stateODEGaussPointHistory}.
""", [("computeStableDeltaT", "El paso estable de referencia.")]),
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


def all_sources():
    """El .C primero, después los headers compartidos: (path, relativo)."""
    out = [(SOLVER / SOURCE, "../" + SOURCE)]
    for h in sorted(COMMON.glob("*.H")):
        out.append((h, "../../../../src/monodomainCommon/" + h.name))
    return out


def main():
    src = SOLVER / SOURCE
    sources = all_sources()
    text = "\n".join(p.read_text() for p, _ in sources)
    located = {}                      # nombre -> (rel, a, b); el .C manda
    per_file = []
    for path, rel_path in sources:
        fx = find_functions(path)
        per_file.append((path.name, fx))
        for name, (a, b) in fx.items():
            located.setdefault(name, (rel_path, a, b))
    funcs = {n: (a, b) for n, (_, a, b) in located.items()}

    body = []
    for title, prose, items in SECTIONS:
        body.append(r"\section{%s}" % title)
        body.append(prose.strip())
        for ref, caption in items:
            if (isinstance(ref, tuple) and isinstance(ref[0], str)
                    and str(ref[1]).startswith("@")):
                # (archivo, "@marcador de inicio", "marcador de fin"): un bloque
                # que no es una función entera, ubicado por texto y no por número
                # de línea, así que el rango se recalcula si el código se mueve.
                fname, start_mark, end_mark = ref
                path, rel = next((pp, rr) for pp, rr in sources if pp.name == fname)
                lines = path.read_text().split("\n")
                try:
                    a = next(i for i, l in enumerate(lines) if start_mark[1:] in l) + 1
                    b = next(i for i in range(a, len(lines)) if end_mark in lines[i])
                except StopIteration:
                    print(f"  aviso: no encontré el bloque {start_mark!r} en {fname}",
                          file=sys.stderr)
                    continue
                trunc = False
            elif isinstance(ref, tuple) and isinstance(ref[0], str):
                # (archivo, nombre, línea desde la que buscar): para nombres que
                # se repiten en un archivo, como los métodos de dos clases.
                fname, name, after = ref
                path, rel = next((pp, rr) for pp, rr in sources if pp.name == fname)
                try:
                    a, b, trunc = span(find_functions(path, after), name)
                except KeyError as exc:
                    print(f"  aviso: {exc}", file=sys.stderr)
                    continue
            elif isinstance(ref, tuple):
                a, b, trunc = ref[0], ref[1], False
                rel = "../" + SOURCE
            else:
                try:
                    a, b, trunc = span(funcs, ref)
                except KeyError as exc:
                    print(f"  aviso: {exc}", file=sys.stderr)
                    continue
                rel = located[ref][0]
            caption = caption + r" \texttt{[%s]}" % rel.split("/")[-1].replace("_", r"\_")
            body.append(emit_listing(rel, a, b, caption, trunc))

    # indice completo de funciones
    index = [r"\section{Índice de funciones}",
             "Todas las funciones, archivo por archivo y en orden de aparición: "
             "primero las del solver, después las de los headers compartidos de "
             r"\texttt{src/monodomainCommon}. Las secciones anteriores recorren "
             "las principales; ésta es el mapa completo."]
    for fname, fx in per_file:
        if not fx:
            continue
        index += [r"\subsection*{\texttt{%s}}" % fname.replace("_", r"\_"),
                  r"\begin{longtable}{rll}", r"\toprule",
                  r"línea & función & líneas \\", r"\midrule", r"\endhead"]
        for name, (a, b) in sorted(fx.items(), key=lambda kv: kv[1][0]):
            if name == "if":
                continue
            index.append(r"%d & \texttt{%s} & %d \\" % (a, name.replace("_", r"\_"), b - a + 1))
        index += [r"\bottomrule", r"\end{longtable}"]

    tex = "\n\n".join([
        PREAMBLE.replace('LITERATE_HERE', literate_clause(text)),
        r"\title{\textbf{Recorrido del código}\\[4pt]\large\texttt{%s}}"
        % NAME.replace("_", r"\_"),
        r"\author{}",
        r"\date{Generado el %s a partir de \texttt{%s} (%d líneas) y de los "
        r"headers de \texttt{src/monodomainCommon}}"
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
