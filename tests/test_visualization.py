"""
Tests for analysis.visualization module
"""
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv

from gdmate import visualization

# Create very small mesh of 9 cells with 16 points and assign each a value
mesh = pv.UniformGrid(dims=(4, 4, 1)).cast_to_unstructured_grid()
mesh['sample_field'] = np.arange(9)

def test_plot2D():
    """ Test plot2D function """
    fig,ax = plt.subplots(1)
    ax = visualization.plot2D(mesh,'sample_field',bounds=[0,2,0,2],ax=ax)

    # Test that the axes contains an image
    assert len(ax.get_images())== 1
    
    # Test that the axes x limits are correct
    assert ax.get_xlim() == (0,2)
    
