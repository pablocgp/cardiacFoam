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


def find_line(lines, needle):
    """Numero de linea (1-based) de la UNICA linea que contiene `needle`.

    Falla si no aparece, y falla tambien si aparece mas de una vez. Lo segundo
    es a proposito: un ancla ambigua elige un sitio en silencio, que es la misma
    falla que un rango escrito a mano -- termina apuntando a otra cosa y nadie
    se entera. Mismo criterio que el chequeo de LITERATE.
    """
    hits = [i + 1 for i, ln in enumerate(lines) if needle in ln]
    if not hits:
        raise KeyError(f"no encontre el ancla {needle!r} en el fuente")
    if len(hits) > 1:
        raise KeyError(f"el ancla {needle!r} aparece {len(hits)} veces "
                       f"(lineas {', '.join(map(str, hits))}); hacela mas "
                       f"especifica")
    return hits[0]


def resolve(ref, lines, funcs, limit=MAX_LISTING):
    """(primera, ultima, recortada?) para una referencia del recorrido.

    `ref` puede ser:
        "nombre"            la funcion, localizada por nombre contando llaves;
        (ancla, ancla)      de una linea de texto a otra, ambas incluidas;
        (ancla, ancla, k)   idem, mas k lineas despues de la segunda;
        (ancla, n)          n lineas a partir del ancla.

    Ya NO se aceptan pares de enteros. Habia tres rangos escritos a mano para
    los fragmentos que no son una funcion entera, y los tres terminaron
    apuntando a otra cosa cuando el fuente se movio -- uno caia en medio del
    comentario de una funcion distinta. Era el mismo problema que copiar el
    codigo dentro del reporte, disfrazado de numero.
    """
    if isinstance(ref, str):
        return span(funcs, ref, limit)

    a = find_line(lines, ref[0])
    if isinstance(ref[1], int):
        b = a + ref[1] - 1
    else:
        b = find_line(lines, ref[1]) + (ref[2] if len(ref) > 2 else 0)

    b = min(b, len(lines))
    if b - a + 1 > limit:
        return a, a + limit - 1, True
    return a, b, False


# ---------------------------------------------------------------------------
# el recorrido
# ---------------------------------------------------------------------------
# (titulo, prosa, [(nombre_funcion | (primera,ultima), pie)])

SECTIONS = [
 ("Qué resuelve este solver", r"""
Resuelve la ecuación de monodominio con una \textbf{solución manufacturada}: se
elige de antemano un campo exacto $V_m(x,t)$ y tres estados $u_1,u_2,u_3$, se los
mete en la ecuación, y lo que sobra se agrega como término fuente. Con eso el
error es \emph{medible} y no estimable, que es lo que permite demostrar el orden
de precisión --- el resultado central del solver.

La ecuación discretizada en el tiempo, con $\theta$ eligiendo el esquema
($\theta=1$ backward Euler, $\theta=1/2$ Crank--Nicolson), es
%
\[ \Big(\tfrac{M}{\Delta t} - \theta\,\lambda K\Big)\,V_m^{n+1}
   \;=\; \Big(\tfrac{M}{\Delta t} + (1-\theta)\,\lambda K\Big)\,V_m^{n}
   \;+\; b \]
%
con $M$ la matriz de masa, $K$ el operador de difusión y $b$ el término fuente
más las condiciones de borde. El archivo está organizado casi exactamente en ese
orden: primero la infraestructura lineal, después $M$ y $K$, después la parte no
lineal, y al final el lazo temporal.

\textbf{Dos ejes de paralelismo}, y conviene tenerlos separados desde el
principio porque el archivo los trata en lugares distintos: \textbf{MPI} reparte
la malla entre ranks y toca todas las fases, mientras \textbf{OpenMP} reparte un
lazo de celdas dentro de un rank y sólo toca la integración de las EDO.
""", []),

 ("De qué se compone: tres archivos y una capa", r"""
Conviene decirlo antes de mirar un listado, porque los nombres que aparecen en
todos ellos no dicen de dónde salen. El solver se construye a partir de
\textbf{tres} archivos del repositorio y de nada más:

\begin{itemize}
  \item \texttt{highOrderManufacturedFDAImplicitPETScDistributed.C}, el recorrido
        de este documento;
  \item \texttt{createFields.H}, que lee los diccionarios y construye los campos
        y los operadores;
  \item \texttt{src/highOrderAdapter/highOrderInterp.H}, la capa de alto orden.
\end{itemize}

Afuera queda OpenFOAM, PETSc, Eigen y \textbf{solids4foam}, de donde sale la
reconstrucción de alto orden propiamente dicha. Después de sacar del
\texttt{Make/options} las inclusiones y bibliotecas que ya no se usaban, el
binario enlaza tres objetos compartidos:
\texttt{libsolids4FoamModels}, \texttt{libpetsc} y \texttt{libOpenFOAM}.

\textbf{Por qué el código dice ``LRE'' en todas partes.} La reconstrucción la
hacía originalmente LRE (Castrillo \emph{et al.}, \emph{Comput Struct}
268:106829, 2022), una biblioteca \textbf{serial}: su construcción de stencils
retorna antes del intercambio MPI --- que está comentado --- y los patches de
procesador lanzan \texttt{NotImplemented}. Es decir que el port a MPI no podía
hacerse sobre ella. La reemplaza \texttt{movingLeastSquares}, también de
solids4foam, que sí es paralela de verdad: los stencils llevan
\textbf{identificadores globales} de celda y cada punto de evaluación hace su
propio intercambio de halo, de modo que el llamador no cambia para correr en
paralelo. Los nombres \texttt{LREInterp\_Vm}, \texttt{LREInterp\_Iion},
\texttt{LREInterp\_states} y las claves \texttt{LRECoeffs*} del diccionario son
\textbf{históricos}: LRE ya no interviene.

\textbf{Qué hace la capa}, porque casi todo lo que contiene es reenvío de una
línea y podría no existir. Cuatro diferencias de contrato la justifican:

\begin{enumerate}
  \item \textbf{Las claves del diccionario.} El solver y todos los casos del
        árbol de corridas están escritos en la ortografía de LRE
        (\texttt{N}, \texttt{Nn}, \texttt{weightFunction}, \texttt{k}), y la capa
        las traduce a las de \texttt{movingLeastSquares}. Vale la pena saber que
        \texttt{N} es el \emph{grado} del polinomio y \texttt{Nn} son las celdas
        \emph{agregadas} sobre el mínimo, que es el número de términos del
        desarrollo de Taylor: $(p{+}1)(p{+}2)/2$ en 2D y
        $(p{+}1)(p{+}2)(p{+}3)/6$ en 3D. Con \texttt{N 2} y \texttt{Nn 10} en 2D
        el stencil tiene $6+10=16$ celdas.
  \item \textbf{Los coeficientes de derivadas por debajo del orden que los
        necesita.} LRE reservaba los arreglos de derivada segunda y tercera
        siempre, así que el solver liga la referencia y recién guarda el
        \emph{uso} con \texttt{order() >= 2}. \texttt{movingLeastSquares} no los
        construye por debajo del orden requerido y \textbf{aborta} si se los
        piden; la capa devuelve una lista vacía y eso es lo que mantiene
        correctos los tres bindings incondicionales del solver.
  \item \textbf{Los pesos de cuadratura cambiaron de significado.} LRE los
        normalizaba --- una cara sumaba 1 --- y el llamador multiplicaba por
        $|S_f|$; los de solids4foam son \textbf{físicos} y suman $|S_f|$. Para
        que la diferencia no reescalara un resultado en silencio, los accesores
        se llaman \texttt{faceQuadWeightPhysical()} y
        \texttt{cellQuadWeightPhysical()}: un sitio portado sin revisar
        \textbf{no compila}.
  \item \textbf{Una sobrecarga rota.} La forma de \texttt{fGrad} que devuelve
        \texttt{autoPtr} llama a un miembro que no existe, así que la capa
        dimensiona el resultado y usa la de dos argumentos.
\end{enumerate}

Nada de eso desaparece si se borra la capa: se muda a cada consumidor. Y hay dos
consumidores, porque el solver del electro usa la misma.
""", []),

 ("Andamiaje de PETSc: errores y matriz distribuida", r"""
El solver ensambla en Eigen y resuelve en PETSc, así que lo primero del archivo
es el puente entre los dos.

\texttt{checkPetscError} convierte un código de retorno de PETSc en un
\texttt{FatalError} de OpenFOAM con el contexto de la llamada. Sin eso un error
de PETSc se propaga como un entero que nadie mira.

\texttt{buildDistributedMat} es la pieza central del port a MPI: arma la matriz
\texttt{MPIAIJ} a partir de los tripletes locales. La matriz es
\textbf{nLocal $\times$ nGlobal} --- filas propias, columnas globales --- que es
el layout que PETSc espera. En serie \texttt{gRowStart} vale cero y
\texttt{nLocal == nGlobal}, así que el camino es el mismo.
""", [("checkPetscError", "Un código de retorno de PETSc pasa a ser un error legible."),
      ("buildDistributedMat", "Tripletes locales a matriz distribuida. Recortado.")]),

 ("El offset global de fila", r"""
\texttt{gRowStart} es el desplazamiento de la primera fila de este rank, o sea
\texttt{globalIndex::localStart()}. Es constante para toda la corrida --- una
malla, una descomposición --- así que se fija una vez en \texttt{main()} en vez
de arrastrarlo por la firma de cada rutina.

Vale la pena detenerse acá porque es donde vivían cuatro bugs del port, todos
del mismo tipo: escribir una columna \emph{local} donde iba una \emph{global}.
En serie \texttt{gRowStart} es cero y \texttt{toGlobal()} es la identidad, así
que \textbf{ninguno de los cuatro era visible sin correr en paralelo}. El test
\textbf{M09} compara $K$ entrada por entrada contra el ensamblado serial
deshaciendo la permutación con \texttt{cellProcAddressing}, y es lo único que
puede verlos: un fingerprint de nnz y sumas de fila es invariante a mover una
entrada de columna.

Justamente por eso el valor no se lee directo sino a través de
\texttt{globalRowStart()}, que \textbf{se niega a servir el cero inicial}. La
asignación ocurre en \texttt{main()}, después de construir el interpolador, y
todas las inserciones de fila en PETSc la leen. Un camino de código que llegara
a una matriz antes de ese punto usaría offset cero \emph{en todos} los ranks, o
sea que los ranks $>0$ escribirían sus filas encima de las del rank 0: operador
corrompido y respuesta plausible pero equivocada, nunca un error. Hoy nada lo
hace; el accesor lo vuelve imposible en vez de meramente improbable.
""", [(("// Global row offset", "return gRowStart;", 1),
       "El offset, y el accesor que se niega a servirlo sin asignar.")]),

 ("Nombres de opciones: una sola ortografía, y los de PETSc", r"""
\textbf{Cada opción del diccionario tiene exactamente una ortografía}, comparada
distinguiendo mayúsculas, y cualquier otra cosa detiene la corrida al arrancar
con la lista de valores válidos. Lo hace \texttt{requireOneOf}, que se llama en
\texttt{createFields.H} justo después de leer cada clave.

Reemplazó a una capa que pasaba todo a minúsculas y acumulaba alias
(\texttt{sparselu}/\texttt{lu}, \texttt{picard}, \texttt{localDiagonal},
\texttt{RKT45}, \texttt{consistentHO}, \texttt{jakobi}\ldots). No era cosmético:
así se escondió una respuesta equivocada. Con \texttt{linearSolver SparseLU}, el
camino de diagonalIion forzaba un solve directo, pero el solver cacheado de
Picard recibía \texttt{preonly} con el PC por defecto, \texttt{ilu}: \emph{una}
aplicación de ILU haciéndose pasar por la inversa. Picard convergía al punto fijo
del operador equivocado --- $V_m$ con un error $4{,}6$ veces mayor, cero pasos
no convergidos, ninguna advertencia. \texttt{resolvePetscKspPc} resuelve ahora
ese par en un solo lugar para los tres caminos PETSc.

Quedan dos vocabularios separados. \texttt{linearSolver} es neutral
(\texttt{GMRES}, \texttt{BiCGSTAB}, \texttt{SparseLU}); las claves de PETSc usan
los nombres propios de PETSc de listas cerradas. \texttt{petscKspTypeForLinearSolver}
es lo único que traduce de uno a otro, para el valor por defecto de
\texttt{petscLinearKspType}.

La función de esta sección que importa para interpretar resultados es
\texttt{petscParallelPcTypeName}.

\texttt{ilu} y \texttt{lu} \textbf{no tienen implementación MPIAIJ}, así que en
paralelo PETSc los envuelve en \textbf{block-Jacobi}, y el bloque \emph{es} la
partición. Es decir que a np=1 y a np=2 no se resuelve con el mismo
precondicionador: son dos algoritmos distintos. \texttt{hypre} y \texttt{gamg}
sí son paralelos nativos, pero su \emph{coarsening} también depende de cómo
quedó cortada la malla.

Sólo \texttt{jacobi} es el mismo operador a cualquier número de ranks, y ésa es
la razón de que la suite lo fije en todos los tests de invariancia. El test
\textbf{M07} lo cuantifica: la desviación entre np=1 y np=2 es
\textbf{1.26e-11} con \texttt{jacobi} contra \textbf{3.60e-05} con \texttt{ilu},
un factor de 2.86 millones.

Para producción la respuesta es otra y está en \textbf{M18}: ahí lo que importa
es el tiempo, y el ranking cambia con el paso de tiempo.
""", [("requireOneOf", "Una sola ortografía por opción, validada al arrancar."),
      ("resolvePetscKspPc", "El par (KSP, PC) de un linearSolver, en un solo lugar."),
      ("petscParallelPcTypeName", "La sustitución que hace que ilu deje de ser ilu en paralelo.")]),

 ("La solución manufacturada", r"""
Acá está lo que hace que este solver pueda demostrar un orden en vez de
estimarlo. \texttt{exactVm} y las tres \texttt{exactU} dan el campo elegido;
\texttt{computeF} y \texttt{computeG} son sus factores espaciales, escritos por
dimensión.

\texttt{vmSourcePDE} es el término que se agrega a la ecuación para que ese campo
sea solución exacta: se obtiene metiendo la solución elegida en el operador y
despejando lo que falta.

\texttt{computeBeta} lee \texttt{conductivity[0]}, o sea la celda cero, y eso
sólo es válido con conductividad \textbf{uniforme}: en una malla descompuesta la
celda cero es una celda física distinta en cada rank. No es un item abierto sino
un \textbf{supuesto asumido} --- este solver existe para este caso manufacturado,
donde $D$ es uniforme por diseño. Queda registrado acá porque es lo primero que
habría que tocar si alguna vez entraran fibras, que hacen $D$ no uniforme por
construcción.

\textbf{Hay dos soluciones manufacturadas}, y la clave
\texttt{manufacturedSolution} elige. La de siempre, \texttt{sines}, es
$\sqrt{1+t}\,\cos(\pi x)\cos(2\pi y)\cdots$: suave y no polinómica, así que
ninguna reconstrucción de acá es exacta sobre ella y el orden \emph{espacial}
medido es el del esquema.

La otra, \texttt{uniform}, tiene todos los factores espaciales iguales a uno.
Existe por una razón puntual: toda reconstrucción, cuadratura y flujo de este
solver es exacto sobre un campo constante, así que el error espacial es cero a
nivel de redondeo y el error contra la solución exacta es \textbf{puramente
temporal}. Es la única forma de ajustar un orden temporal contra la solución
analítica en vez de contra otra solución numérica, y es la vía B de \textbf{M02}.
A malla fija con \texttt{sines} eso es imposible: el error espacial
($2{,}4\times10^{-5}$ en hexa $N=40$) entierra al temporal, que para BDF4 es
$\sim10^{-10}$.

Lo que cuesta, dicho para que nadie lea de más el resultado: el laplaciano de una
constante es cero, así que bajo \texttt{uniform} la difusión no participa y lo
que se verifica es la integración temporal del sistema de reacción. La difusión
la cubre \textbf{M02} por autoconvergencia. Y \texttt{computeBeta} devuelve cero
en ese modo, porque no hay laplaciano que cancelar.
""", [("exactVm", "El campo exacto."),
      ("computeF", "El factor espacial, y dónde el modo uniform lo vuelve 1."),
      ("vmSourcePDE", "El término que lo convierte en solución."),
      ("computeBeta", "Lee la celda 0: sólo vale con D uniforme.")]),

 ("Bloques de construcción del ensamblado", r"""
Antes de las dos matrices de rigidez vienen las piezas que las dos comparten.

\texttt{addTripletIfNeeded} es el punto por el que pasa \emph{toda} entrada de la
matriz. Descarta las que valen menos que \texttt{SMALL}, lo que hace que el
patrón de esparsidad dependa del valor y, por lo tanto, de la partición al nivel
del redondeo.

\texttt{addCellGradientDotCoeffs} agrega la fila de gradiente de una celda
proyectada sobre un vector. Usa \texttt{globalCellStencils()}, y ahí vivió durante
meses el que se creía el problema abierto más grande del proyecto: sobre las
celdas pegadas a un corte de procesador el stencil \textbf{perdía un miembro} (17
en serie contra 16 en paralelo). Como el ajuste MLS se rehace sobre el conjunto
reducido, no cambiaba sólo el término que faltaba --- \textbf{refitaba toda la
fila}. Eso era lo que hacía fallar a \textbf{M11} y, río abajo, a \textbf{M15}, y
se había dado por un item de solids4foam que se mide y no se arregla.

\textbf{No lo era.} La causa estaba del lado del caso: la clave \texttt{Nn} del
diccionario son las celdas que se AGREGAN sobre el mínimo que exige el grado, y
los barridos la venían escribiendo como si fuera el TOTAL. La biblioteca sumaba
el mínimo otra vez, así que los estencils medían $2\,\mathrm{min}+10$ en vez de
$\mathrm{min}+10$ --- 16/22/30 celdas en 2D en vez de 13/16/20, un 50 \% más
grandes de lo pretendido --- y el halo no alcanzaba a cubrirlos. Con el tamaño
correcto el stencil es invariante a la partición: \textbf{M11 mide 0 stencils
distintos y M15 baja de 4.981e-06 a 6.953e-15}, los dos en PASS. La suite quedó
24/24.

\texttt{exchangeCoupledFaceRows} es la que permite que el término de
estabilización cruce un corte: intercambia con el rank vecino las filas que
hacen falta del otro lado. Es \textbf{colectiva}, y está deliberadamente fuera
del lazo de chunks --- adentro se trabaría en cuanto dos ranks tuvieran una
cantidad distinta de chunks.
""", [("addTripletIfNeeded", "Por acá pasa toda entrada de la matriz."),
      ("addCellGradientDotCoeffs", "La fila de gradiente. Recortada."),
      ("exchangeCoupledFaceRows", "Intercambio de filas a través de un corte. Recortada.")]),

 ("Las dos matrices de rigidez", r"""
El solver tiene dos ensamblados del operador de difusión y la tríada del caso
elige cuál corre.

\texttt{assembleStandardOrthogonalStiffnessMatrix} es el de bajo orden: flujo de
dos puntos entre celdas vecinas. Usa el \textbf{promedio}
$\tfrac12(D_{own}+D_{nei})$ como tensor de cara.

\texttt{assembleHighOrderStiffnessMatrix} es el de alto orden, sobre el stencil
MLS con cuadratura en las caras. Acá estuvo el defecto más grande que encontró
la verificación: usaba \texttt{conductivity[own]} \textbf{en las dos filas} de
una cara interna, en vez del promedio. Con $D$ uniforme los dos tensores
coinciden y se cancela idénticamente, así que fue invisible durante toda la vida
del código; con $D$ no uniforme hacía que $K$ \textbf{cambiara un 21\,\% al
descomponer la malla}, porque al cortar, cada rank es dueño de su lado de la
cara. Está corregido, y \textbf{M15} lo mide: de 2.108e-01 a 4.736e-15.

Las dos ensamblan \textbf{por chunks} de caras, y la clave del diccionario
\texttt{assemblyFaceChunk} controla el tamaño. Es control de capacidad
solamente: el resultado no debe depender de él, y \textbf{M12} lo verifica ---
con tolerancia y no bit a bit, porque cambiar el tamaño de bloque cambia el
orden de suma en punto flotante.
""", [("assembleStandardOrthogonalStiffnessMatrix", "Bajo orden: promedio en la cara. Recortado."),
      ("assembleHighOrderStiffnessMatrix", "Alto orden, sobre el stencil MLS. Recortado.")]),

 ("Resolución del sistema lineal", r"""
\texttt{solveSparseSystem} despacha entre PETSc y Eigen según
\texttt{linearSolverBackend}. En paralelo sólo el camino de PETSc sirve, y ahora
eso \textbf{se hace cumplir} en vez de sólo estar comentado: la matriz es
\texttt{nLocal $\times$ nGlobal}, o sea rectangular, y las factorizaciones de
Eigen necesitan una matriz cuadrada. Sin el guard, \texttt{SparseLU} moría en un
\texttt{eigen\_assert} que nombra un archivo interno de la librería --- ni este
solver, ni la clave del diccionario que lo eligió, ni la palabra ``paralelo''.
Medido a np=2: el error ahora dice \texttt{200 local rows and 400 global
columns}.

Hay además un GMRES escrito a mano (\texttt{solveGMRES}) que queda como camino
de respaldo cuando no se usa PETSc. Ése tiene su propio guard, y por un motivo
\textbf{distinto}: no hay matriz que pueda ser rectangular, pero sus productos
internos (\texttt{V[i].dot(w)}, \texttt{w.norm()}) los calcula Eigen sobre el
vector \emph{local} y nunca se reducen. Cada rank armaría su propia base de
Krylov y su propio test de convergencia, mientras el producto matriz-vector que
llama \emph{sí} es colectivo (\texttt{MatMult} de PETSc). Dos ranks que salgan
del lazo en distinta iteración se traban ahí: no es un abort, es un cuelgue.

Un tercer punto de la misma familia: \texttt{usesPetscBackend} \textbf{valida} el
nombre en vez de tratarlo como una lista blanca de dos entradas cuyo
\texttt{else} era Eigen. Cualquier valor que no fuera \texttt{petsc}/\texttt{ksp}
--- un typo en un diccionario, por ejemplo --- mudaba la corrida entera al
backend serial sin una línea en el log. Es la misma falla que
\texttt{nonlinearMethod Picrad} en el solver electro, que colapsaba una matriz de
métodos sobre un solo método.

\textbf{Qué solver de Krylov corresponde} lo decide la simetría del operador, y
está medida en \textbf{M10}: el de bajo orden es simétrico a 8.4e-16, el de alto
orden \textbf{no} --- 7.9\,\% de asimetría relativa, con 300 de 4292 entradas sin
transpuesta. Con la tríada de alto orden hay que usar \texttt{gmres}.
\texttt{cg} no falla siempre ahí, que es lo que lo hace peligroso: sobrevive
mientras le alcancen pocas iteraciones y se cae cuando necesita más
(\textbf{E13} lo muestra fallando con \texttt{jacobi} y \texttt{gamg}).
""", [("usesPetscBackend", "El nombre del backend, validado."),
      ("solveSparseSystem", "El despacho entre backends, con el guard de paralelo. Recortado."),
      ("solveSparseSystemEigen", "El camino serial. Recortado.")]),

 ("Discretización temporal: las dos familias", r"""
El solver ofrece dos familias de integrador temporal, y \texttt{implicitScheme}
elige entre ellas. Con el orden que cada una debe entregar:

\begin{center}
\begin{tabular}{@{}lll@{}}
\toprule
\texttt{implicitScheme} & orden & estabilidad \\
\midrule
\texttt{backwardEuler} & $O(\Delta t)$ & A y L-estable \\
\texttt{crankNicolson} & $O(\Delta t^2)$ & A-estable, \textbf{no} L-estable \\
\texttt{BDF1} & idéntico a \texttt{backwardEuler}, bit a bit & A y L \\
\texttt{BDF2} & $O(\Delta t^2)$ & A y L-estable \\
\texttt{BDF3} & $O(\Delta t^3)$ & $A(86{,}03^\circ)$ \\
\texttt{BDF4} & $O(\Delta t^4)$ & $A(73{,}35^\circ)$ \\
\bottomrule
\end{tabular}
\end{center}

\texttt{timeSchemeCoeffs} traduce el nombre a los coeficientes. Lo que hay que
ver en esa función es que \textbf{las dos familias terminan en el mismo par de
operadores}: $A=(a_0/\Delta t)M-\theta L$ y $B=M/\Delta t+(1-\theta)L$. BDF
cambia \emph{sólo} el escalar que multiplica a $M$, así que el ensamblado, el
patrón de esparsidad y la factorización cacheada quedan intactos. Y el historial
entra por un \emph{único} producto por la matriz de masa, porque la combinación
$H^n=-\sum_{j\ge1}a_jV^{n+1-j}$ se arma como campo antes de aplicar $B$: un paso
BDF-$k$ cuesta un producto matriz--vector, no $k$.

BDF5 y BDF6 no se ofrecen. Pasada la segunda barrera de Dahlquist ningún
multipaso es A-estable, y la cuña se cierra de $86^\circ$ y $73^\circ$ en BDF3/4
a $51{,}8^\circ$ y $17{,}8^\circ$ en BDF5/6.

\textbf{El arranque es la parte que engaña.} Un método de $k$ pasos necesita
valores iniciales precisos a $O(\Delta t^k)$ \emph{respecto del sistema
semi-discreto}, no de la PDE. Sembrar con la solución manufacturada da
$O(\tau_h\,\Delta t)$ --- primer orden --- porque la solución exacta de la PDE
no satisface las ecuaciones semi-discretas. Nada avisa: la corrida termina bien y
sólo el orden sale mal, que es por lo que \texttt{bdfStartup} tiene hoy
\texttt{ESDIRK3} \textbf{por defecto}.

Con ese arranque los órdenes son los nominales. Medido a $N=40$ por las dos vías
que usa la suite --- \textbf{M02}, autoconvergencia de Cauchy sobre la solución
de senos con la difusión activa, y su vía B, error contra la solución
analítica sobre la MMS uniforme:

\begin{center}
\begin{tabular}{@{}llll@{}}
\toprule
\texttt{implicitScheme} & nominal & Cauchy, vía A & analítica, vía B \\
\midrule
\texttt{backwardEuler} y \texttt{BDF1} & 1 & $1{,}000$ & $1{,}000$ \\
\texttt{crankNicolson} & 2 & $1{,}994$ & $2{,}000$ \\
\texttt{BDF2} & 2 & $1{,}979$ & $1{,}992$ \\
\texttt{BDF3} & 3 & $2{,}932$ & $2{,}971$ \\
\texttt{BDF4} & 4 & $3{,}890$ & $3{,}941$ \\
\bottomrule
\end{tabular}
\end{center}

\texttt{exact} y \texttt{constant} se conservan, pero \textbf{sólo} para que el
costo de un mal arranque siga siendo medible: es lo que hace \textbf{M20},
aislando el efecto como
$e(\Delta t)=\|V_{\text{arranque}}-V_{\text{ESDIRK3}}\|$ al mismo $\Delta t$ y la
misma malla, donde el error propio de BDF3 se cancela. Los dos modos resultan
$O(\Delta t)$, y difieren donde la teoría dice: el de \texttt{constant} no cambia
al refinar la malla (razón $1{,}00$) y el de \texttt{exact} cae $15{,}5$ veces
por refinamiento, porque su coeficiente es $\tau_h$.

Por eso existe \texttt{esdirk3Tableau}, y \textbf{sólo} por eso: ESDIRK no es
seleccionable como esquema. Es autoarrancable, que es lo que el arranque pide,
pero tiene orden de etapa 2 y sufre reducción de orden --- ESDIRK4 midió $3{,}03$
donde BDF4 llega a $3{,}89$; esa cifra de ESDIRK4 es anterior al arreglo del lazo
temporal y no se volvió a medir --- y cada etapa reintegra las EDO, que es el
$74$--$99\%$ del tiempo en el solver fisiológico.
""", [("timeSchemeCoeffs", "Los coeficientes de las dos familias. Recortado."),
      ("esdirk3Tableau", "La tabla del arranque, derivada de las condiciones de orden."),
      (("// The ESDIRK startup step",
        "const bool esdirkThisStep =", 1),
       "Dónde el arranque ESDIRK es alcanzable, y sólo ahí.")]),

 ("El driver de Vm que manejan las EDO", r"""
Las EDO de estado se integran con sub-pasos \emph{dentro} del paso de la PDE, y
en cada sub-paso necesitan un $V_m$ en un instante intermedio que el esquema no
conoce. Cómo se rellena ese hueco \textbf{acota el orden de todo el solver}: una
recta entre los dos extremos es $O(\Delta t^2)$, y ningún integrador temporal por
encima de segundo orden puede mostrar el suyo detrás de eso.

No es una salvedad teórica, está medido. La clave \texttt{stateODEDriver} conmuta
entre las dos opciones --- \texttt{history}, el interpolante por los niveles que
el esquema ya guarda, y \texttt{linear}, la rampa vieja --- y \textbf{M21} corre
BDF4 con cada una sin cambiar nada más: la rampa ajusta $\mathbf{2{,}04}$ contra
$\mathbf{3{,}89}$ del interpolante, con el error del par más fino mayor por un
factor $44{,}7$. El $2{,}04$ es exactamente lo que predice la cuenta: la rampa le
erra por $O(\Delta t^2)$ punto a punto, que integrado sobre el paso es
$O(\Delta t^3)$ local y $O(\Delta t^2)$ global.

\texttt{history} es el valor por defecto; \texttt{linear} existe para que ese
costo siga siendo medible, porque un orden capado es todo lo que la corrida
reporta --- no hay error ni advertencia. El camino \texttt{gaussPointODE} usa
siempre la rampa: subirlo obligaría a guardar cada nivel \emph{en cada punto de
cuadratura}.

\texttt{interpolateVmDriver} evalúa el polinomio de Lagrange a través de los
voltajes que el esquema exterior ya calculó: dos nodos para las familias
$\theta$ y BDF1, $k+1$ niveles de historial para BDF-$k$, y los valores de etapa
para el arranque ESDIRK. No cuesta ninguna evaluación extra --- los nodos ya
existen. El caso de dos nodos está escrito aparte y no pasa por el bucle, que es
lo que garantiza que \texttt{backwardEuler} y \texttt{crankNicolson} sigan
dando resultados idénticos bit a bit.
""", [("interpolateVmDriver", "El interpolante, con el caso de dos nodos escrito aparte."),
      ("reactionRatesLinearVm", "Dónde el driver entra a las EDO.")]),

 ("Reconstrucción a los puntos de integración", r"""
Con alto orden, $V_m$ y los estados viven en los centros de celda pero la
corriente iónica se evalúa en los puntos de Gauss. Estas rutinas hacen la
reconstrucción MLS de ida y el promedio de vuelta.

\texttt{reconstructStatesAtIionIntegrationPoints} lleva un limitador de
Barth--Jespersen, acotando la reconstrucción al rango de los vecinos. No es
cosmético: sin él, el sobrepaso de alto orden sobre estados rígidos desborda.
""", [("reconstructVmAtIionIntegrationPoints", "Reconstrucción de Vm. Recortada."),
      ("averageIntegrationPointFieldToCells", "El promedio de vuelta a celdas.")]),

 ("Integración de las EDO, y el paralelismo por hilos", r"""
\texttt{advanceStateODE} avanza los estados en cada celda o punto de integración,
con RK4 o RKF45 adaptativo según el diccionario.

Acá están los \textbf{únicos} lazos con OpenMP de todo el archivo. Es importante
para entender el rendimiento: los hilos reparten \emph{estas} iteraciones y nada
más --- \textbf{el ensamblado no está paralelizado con hilos}, sólo por MPI. La
región funciona porque cada punto integra su propia EDO sin hablar con ningún
otro, y por eso el resultado es \textbf{bit-idéntico} con cualquier cantidad de
hilos, cosa que \textbf{E20} mide en el otro solver dando exactamente cero.

El umbral por debajo del cual la región no se abre llega como parámetro desde la
clave \texttt{stateODEOpenMPThreshold}, en las \textbf{22} llamadas que lo pasan.
Esto figuró un tiempo como item abierto --- se creía que el valor estaba fijo en
256 en seis lazos y que la clave no hacía nada --- y al ir a arreglarlo resultó
que ya estaba conectado y que lo desactualizado era el registro. Lo que sí se
hizo fue \textbf{quitar los argumentos por defecto} de las seis firmas que
reciben el umbral, para que una llamada nueva que se olvide de pasarlo no
compile, en vez de correr en silencio con un 256 escrito en la firma.
""", [("advanceStateODE", "El avance de los estados. Recortado."),
      ("rkf45StateStep", "El paso adaptativo. Recortado.")]),

 ("main(): montaje", r"""
\texttt{main()} arranca leyendo la malla y los diccionarios, construye los
operadores de alto orden, y ensambla $M$ y $K$ \textbf{una sola vez} --- son constantes,
porque la no linealidad está en el término iónico y no en la difusión. Por eso
en el desglose de tiempos el \texttt{setup} aparece separado del \texttt{loop}.

Es también donde se fija \texttt{gRowStart}, y donde viven los dos volcados de
diagnóstico que la verificación usa, los dos detrás de una variable de entorno y
gratis cuando están apagados: \texttt{CF\_DUMP\_K} escribe $K$ como tripletes
(fila global, columna global, valor) y \texttt{CF\_DUMP\_STENCILS} escribe la
membresía de los stencils. Son la base de M09, M10, M11 y M15.

Lo primero de \texttt{main()}, antes de cualquier asignación grande, son tres
llamadas a \texttt{mallopt}, y no son cosmética: la construcción de los
operadores corre unos \textbf{2.3 millones de factorizaciones QR} sobre una malla
de tetraedros 3D con $N{=}40$ y p3, y sin ajustar el asignador glibc crece el
heap en varios GB de fragmentos que no devuelve y la corrida muere por falta de
memoria en una máquina de 16 GB. \texttt{M\_ARENA\_MAX=2} evita que glibc cree
una arena por hilo de OpenMP --- el montaje es serial, así que las arenas sólo
inflan el RSS ---, y los dos umbrales en 64\,KB hacen que las asignaciones
medianas pasen por \texttt{mmap} y que el relleno del heap vuelva al sistema
operativo. Es el mismo techo de memoria que \textbf{M17} mide y el que explica
por qué en 3D no se pasa de $N{=}30$.
""", [(("int main(int argc, char* argv[])", 24), "El arranque."),
      (("Tighten the glibc allocator before any large allocation. The LRE", 15),
       "El asignador, ajustado antes de la primera asignación grande."),
      (("gRowStart = LREInterp_Vm.globalCells().localStart();", 8),
       "Dónde se fija el offset global de fila."),
      ("logMemoryCheckpoint", "Los checkpoints de memoria que M17 usa para ajustar el techo.")]),

 ("main(): el lazo temporal y los métodos no lineales", r"""
Cada paso resuelve un sistema no lineal, y hay tres métodos: \textbf{Picard} con
relajación, \textbf{JFNK} --- Newton sin jacobiano explícito, con el producto
matriz-vector aproximado por diferencias --- y \textbf{diagonalIion}, que
linealiza la corriente iónica por su derivada diagonal. \textbf{M08} verifica
que los tres lleguen al mismo lugar, y lo que exige no es una cota absoluta sino
una \emph{dirección de cambio}: que la dispersión entre los tres caiga en cada
escalón de una escalera de \texttt{implicitTolerance}. Eso dice exactamente lo
que se quiere afirmar --- que lo que los separa es el piso del solve lineal y no
una diferencia estructural --- y no lleva ninguna constante. Su versión anterior
derivaba una cota de \texttt{nonlinearVmTolerance}, y estaba indexada a una
perilla que no gobernaba la cantidad que acotaba: apretarla cuatro órdenes no
movía la discrepancia ni en el décimo dígito.

Ese test tenía además un punto ciego, y detrás había un defecto real en JFNK. Su
lazo evalúa la convergencia dos veces por iteración, al tope y al final, y la del
tope necesita los incrementos que acumuló la iteración \emph{anterior}. Los
recalculaba en el lugar contra una copia tomada unas líneas antes, con una sola
re-evaluación del \emph{mismo} punto en medio, así que valían cero por
construcción: cuatro de las cinco tolerancias comparaban cero contra un número
positivo y el test se reducía al residual acoplado. Consecuencia medida en 3-D:
en 24 de 50 pasos JFNK devolvía la extrapolación inicial sin dar un paso de
Newton y con cero llamadas al KSP. Se arregló guardando los incrementos medidos
al final de cada iteración; \texttt{max coupled\_relL2} cayó de 9.493e-09 a
5.655e-16 y Picard y \texttt{diagonalIion} quedaron bit a bit idénticos.

Al agotar las iteraciones sin converger hay dos salidas posibles, y \textbf{los
tres métodos las eligen ahora con la misma clave}: \texttt{false} descarta el
paso y vuelve a los valores del inicio, \texttt{true} se queda con el iterado y
lo dice. Antes Picard tenía su propia rama que aceptaba el último iterado
\emph{siempre}, así que el mismo solver tenía dos semánticas de error según el
método. Se eliminó esa rama en vez de agregarle un rollback propio: las otras
dos ya implementan las dos salidas, así que Picard cae solo en la correcta.

Medido forzando la no convergencia con \texttt{nonlinearIterations 1} y
tolerancias en 1e-30, que es la única forma de hacerlo disparar --- en los
barridos de convergencia no ocurre nunca: con \texttt{false} el error final es
\textbf{4.998e-04}, que es el de no haber avanzado el paso, y con \texttt{true}
es \textbf{5.023e-06}. El rollback efectivamente revierte.

Un nombre de método inválido es \texttt{FatalError} y no un silencioso repliegue
sobre Picard, que es lo que hacía el solver electro hasta que \textbf{E18} lo
señaló.

\textbf{La condición de parada del lazo tenía un error, y era invisible.} El
tiempo acumulado es una suma en punto flotante de $\Delta t$, y puede caer por
debajo de \texttt{endTime} \emph{por más que} \texttt{SMALL}: 160 pasos de
$3{,}2\times10^{-3}$ suman $0{,}51199999999999868$, a $1{,}3\times10^{-15}$ de
$0{,}512$. Con la condición vieja, \texttt{endTime - SMALL}, el lazo daba un paso
de más y la corrida terminaba en $0{,}5152$ sin avisar nada. Eso no es cosmético
para un estudio temporal: los campos que se comparan entre corridas con distinto
$\Delta t$ dejan de ser del mismo instante. Ahora la condición es
\texttt{endTime} $-\ \Delta t/2$, exacta porque acá $\Delta t$ es constante. Lo
encontró la guarda de un test de la suite, no una revisión del código, y el mismo
error sigue presente en el solver electro.
""", [(("    const bool usePicard =",
        "Valid options are Picard, JFNK and diagonalIion (exact", 3),
       "El despacho de método no lineal, y el rechazo de un nombre inválido."),
      (("// Picard used to have",
        "nonlinearRolledBack = true;", 1),
       "Las dos salidas al agotar las iteraciones, ahora comunes a los tres métodos."),
      (("// Stop within half a step of endTime",
        "while (runTime.value() <", 1),
       "El lazo que ya no da un paso de más, y por qué."),
      (("// The increments measured at the END of the previous corr",
        "scalar IionIncrPrevIter = GREAT;", 1),
       "Por qué los incrementos del chequeo del tope de JFNK viven fuera del "
       "lazo: recalculados adentro valían cero y mataban cuatro tolerancias.")]),
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
    lines = text.split("\n")
    funcs = find_functions(src)
    rel = "../" + SOURCE

    body = []
    for title, prose, items in SECTIONS:
        body.append(r"\section{%s}" % title)
        body.append(prose.strip())
        for ref, caption in items:
            try:
                a, b, trunc = resolve(ref, lines, funcs)
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
