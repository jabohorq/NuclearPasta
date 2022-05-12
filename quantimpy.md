# QuantImPy
## Python library for image procesing

Using this library, we can calculate the Minkowski functional for MD configurations.
[Documentation](https://boeleman.github.io/quantimpy/source/quantimpy.minkowski.html)

1. Install

`pip install quantimpy`

2. Calculate Minkowski functionals (3D)

```python
# load library
from quantimpy import minkowski as mk
...
res = np.array([2, 2, 2])
minkowski = mk.functionals(image, res)
```

Here `image` is a 3D NumPy boolean array, for example, for a sphere of radius 48 pixels:
```python
from skimage.morphology import (ball)
image = np.zeros([128,128,128],dtype=bool)
image[16:113,16:113,16:113] = ball(48,dtype=bool)
```

Also, `res` is an 3 elements array, with the resolution on each direction. Resolution is defined as 
<unit of length>/<pixel> (lenght per pixel). By default the value is 1 on each direction.
In the example the resolution is 2 on each direction, this means: 2 length units per pixel.
