PyArcKnot
===

PyArcKnot is a package for studying the arc diagrams of mathematical knots.

It was developed as part of my Final Project for my Industrial Design Engineering Degree, with the help and support of Pedro González Manchón.
The Final Project paper (in Spanish) can be found [here](https://oa.upm.es/77063/).


The following is just a showcase of some of the most important functions.
For a more in-depth demonstration of the possibilities available, check out the [JupyterNotebook demo].

The package consists of 4 modules:

```python
import pyarcknot.knot_matrix as km
import pyarcknot.knot_diagram as kd
import pyarcknot.knot_calculate as kc
import pyarcknot.turaev_surface as ts
```

### Knot Matrix

Knot Matrix is used to define the knot diagrams as [NumPy](https://numpy.org/) arrays in order to be able to work with them

```python
k8_21_arc = km.clean_k_arc('2 7 1 4 3 5 4 8 2 6 1 5 3 7 6 8')
k8_21_xco = km.xco_arc(k8_21_arc)
print(k8_21_xco)
```
	[[0 0 0 0 0 1 0 2]
	 [0 0 1 0 0 3 2 0]
	 [1 0 3 0 2 0 0 0]
	 [0 1 3 0 3 2 0 0]
	 [0 0 0 2 3 0 3 1]
	 [0 0 2 3 1 0 0 0]
	 [2 3 0 1 0 0 0 0]
	 [0 2 0 0 0 0 1 0]]


### Knot Diagram

Knot Diagram uses [Matplotlib](https://matplotlib.org/) to display diagrams of the knot.

```python
kd.draw_arc(k8_21_xco)
```

![arc_diagram](https://github.com/inigorrix/pyarcknot/blob/main/docs/arc_diagram.png?raw=true)


It can also return the number of loops in a smoothed diagram

```python
kd.draw_diagrams(k8_21_xco)
```

    Number of crossings = 8

![arc_diagram](https://github.com/inigorrix/pyarcknot/blob/main/docs/smooth_a.png?raw=true)

    |s_A D| = 5 

![arc_diagram](https://github.com/inigorrix/pyarcknot/blob/main/docs/arc_diagram.png?raw=true)

![arc_diagram](https://github.com/inigorrix/pyarcknot/blob/main/docs/smooth_b.png?raw=true)

    |s_B D| = 1 


### Knot Calculate

Knot Calculate is used to calculate properties of the diagram such as the Kauffman Bracket Polynomial or the Jones Polynomial.
For this, [Sympy](https://www.sympy.org/) was used to work with polynomials and [Numba](https://numba.pydata.org/) to increase the performance and reduce the calculation time.

```python
kc.kauffman_bracket(k8_21_xco)
```

$\displaystyle A^{16} - 2 A^{12} + 2 A^{8} - 3 A^{4} + 3 - \frac{2}{A^{4}} + \frac{2}{A^{8}}$

```python
kc.jones_polynomial(k8_21_xco)
```

$\displaystyle \frac{2}{t} - \frac{2}{t^{2}} + \frac{3}{t^{3}} - \frac{3}{t^{4}} + \frac{2}{t^{5}} - \frac{2}{t^{6}} + \frac{1}{t^{7}}$


### Turaev Surface

Turaev Surface uses [Matplotlib](https://matplotlib.org/) to render an image of the a 3D surface obtained from the diagram of the knot.

```python
ts.turaev_surf(k8_21_xco)
```

![arc_diagram](https://github.com/inigorrix/pyarcknot/blob/main/docs/turaev_surface.png?raw=true)
