import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from mayavi import mlab



def plot_unmatched_cycle(BM, input, dim, index, plot_critical_values=False, use='plotly'):
    vertices = BM.get_unmatched_cycle(input, dim, index)

    if (len(BM.shape) == 3):
        cycle = np.zeros(shape=(BM.shape[0]+2, BM.shape[1]+2, BM.shape[2]+2))
        for vertex in vertices[0:-1]:
            cycle[vertex[0]+1, vertex[1]+1, vertex[2]+1] = 1
        if (plot_critical_values):
            if (len(vertices) != 0):
                if (use == 'plotly'):
                    cycle[vertices[0][0]+1, vertices[0][1]+1, vertices[0][2]+1] = 1.1
                    cycle[vertices[-1][0]+1, vertices[-1][1]+1, vertices[-1][2]+1] = 0.9
                if (use == 'mayavi'):
                    cycle[vertices[0][0]+1, vertices[0][1]+1, vertices[0][2]+1] = 2
                    cycle[vertices[-1][0]+1, vertices[-1][1]+1, vertices[-1][2]+1] = -1
        
        if (use == 'plotly'):
            x = np.linspace(1, BM.shape[0]+1, BM.shape[0]+2)
            y = np.linspace(1, BM.shape[1]+1, BM.shape[1]+2)
            z = np.linspace(1, BM.shape[2]+1, BM.shape[2]+2)
            X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

            fig = go.Figure(data=go.Volume( 
                x=X.flatten(), 
                y=Y.flatten(), 
                z=Z.flatten(), 
                value=cycle.flatten(),
                isomin=0.9,
                isomax=1.1,
                opacity=0.1,
                surface_count=100,
                caps= dict(x_show=True, y_show=True, z_show=True),
                ))
            fig.show()

        if (use == 'mayavi'):
            fig = mlab.figure(size=(1000, 1000))
            contour = mlab.contour3d(cycle, colormap='bone', transparent=True)
            mlab.xlabel('X')
            mlab.ylabel('Y')
            mlab.zlabel('Z')
            mlab.show()

    return


def plot_matched_cycles(BM, dim, index, plot_critical_values=False, use='plotly'):
    vertices0, vertices1 = BM.get_matched_cycles(dim, index)

    if (len(BM.shape) == 3):
        cycle0 = np.zeros(shape=(BM.shape[0]+2, BM.shape[1]+2, BM.shape[2]+2))
        for vertex in vertices0[0:-1]:
            cycle0[vertex[0]+1, vertex[1]+1, vertex[2]+1] = 1
        if (plot_critical_values):
            if (len(vertices0) != 0):
                if (use == 'plotly'):
                    cycle0[vertices0[0][0]+1, vertices0[0][1]+1, vertices0[0][2]+1] = 1.1
                    cycle0[vertices0[-1][0]+1, vertices0[-1][1]+1, vertices0[-1][2]+1] = 0.9
                if (use == 'mayavi'):
                    cycle0[vertices0[0][0]+1, vertices0[0][1]+1, vertices0[0][2]+1] = 2
                    cycle0[vertices0[-1][0]+1, vertices0[-1][1]+1, vertices0[-1][2]+1] = -1
        cycle1 = np.zeros(shape=(BM.shape[0]+2, BM.shape[1]+2, BM.shape[2]+2))
        for vertex in vertices1[0:-1]:
            cycle1[vertex[0]+1, vertex[1]+1, vertex[2]+1] = 1
        if (plot_critical_values):
            if (len(vertices1) != 0):
                if (use == 'plotly'):
                    cycle1[vertices1[0][0]+1, vertices1[0][1]+1, vertices1[0][2]+1] = 1.1
                    cycle1[vertices1[-1][0]+1, vertices1[-1][1]+1, vertices1[-1][2]+1] = 0.9
                if (use == 'mayavi'):
                    cycle1[vertices1[0][0]+1, vertices1[0][1]+1, vertices1[0][2]+1] = 2
                    cycle1[vertices1[-1][0]+1, vertices1[-1][1]+1, vertices1[-1][2]+1] = -1
        
        if (use == 'plotly'):
            x = np.linspace(0, BM.shape[0]+1, BM.shape[0]+2)
            y = np.linspace(0, BM.shape[1]+1, BM.shape[1]+2)
            z = np.linspace(0, BM.shape[2]+1, BM.shape[2]+2)
            X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

            fig0 = go.Figure(data=go.Volume( 
                x=X.flatten(), 
                y=Y.flatten(), 
                z=Z.flatten(), 
                value=cycle0.flatten(),
                isomin=0.9,
                isomax=1.1,
                opacity=0.1,
                surface_count=100,
                caps= dict(x_show=True, y_show=True, z_show=True),
                ))
            fig1 = go.Figure(data=go.Volume( 
                x=X.flatten(), 
                y=Y.flatten(), 
                z=Z.flatten(), 
                value=cycle1.flatten(),
                isomin=0.9,
                isomax=1.1,
                opacity=0.1,
                surface_count=100,
                caps= dict(x_show=True, y_show=True, z_show=True),
                ))
            fig = make_subplots(
            rows=1, cols=2,
            specs=[[{'type': 'volume'}, {'type': 'volume'}]])
            fig.add_trace(fig0.data[0], row=1, col=1)
            fig.add_trace(fig1.data[0], row=1, col=2)
            fig.show()

        if (use == 'mayavi'):
            cycle = np.concatenate([cycle0, cycle1], axis=0)
            fig = mlab.figure(size=(1000, 1000))
            contour = mlab.contour3d(cycle, colormap='bone', transparent=True)
            mlab.xlabel('X')
            mlab.ylabel('Y')
            mlab.zlabel('Z')
            mlab.show()

    return