import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from mayavi import mlab



def plot_image(image):
    if len(image.shape) == 2:
        plt.figure(figsize=(4,4))
        plt.imshow(image, cmap='gray')
        plt.axis('off')

    if len(image.shape) == 3:
        fig = mlab.figure(size=(1000, 1000))
        contour = mlab.contour3d(image, colormap='bone', transparent=True)
        mlab.xlabel('X')
        mlab.ylabel('Y')
        mlab.zlabel('Z')
        mlab.show()
    
    return



def plot_unmatched_cycle(BM, input, dim, index, plot_critical_values=False, use='mayavi'):
    vertices = BM.get_unmatched_cycle(input, dim, index)

    if (len(vertices) == 0):
        return
    
    if (len(BM.shape) == 2):
        cycle = np.zeros(shape=(BM.shape[0], BM.shape[1]))
        for vertex in vertices[0:-1]:
            cycle[vertex[0], vertex[1]] = 1
        if (plot_critical_values):
            if (len(vertices) != 0):
                cycle[vertices[0][0], vertices[0][1]] = 1.5
                cycle[vertices[-1][0], vertices[-1][1]] = 0.5
        
        fig = plt.figure(figsize=(15,5))
        fig.add_subplot(1, 1, 1)
        if (plot_critical_values):
            plt.imshow(cycle, cmap='gray', vmin=0, vmax=1.5)
        else:
            plt.imshow(cycle, cmap='gray', vmin=0, vmax=1)
        plt.axis('off')

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

        if (use == 'mayavi'):
            fig = mlab.figure(size=(1000, 1000))
            contour = mlab.contour3d(cycle, colormap='bone', transparent=True)
            mlab.xlabel('X')
            mlab.ylabel('Y')
            mlab.zlabel('Z')
            mlab.show()

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
                surface_count=20,
                caps= dict(x_show=True, y_show=True, z_show=True),
                ))
            #fig.show()
            fig.write_html("plot.html")

    return



def plot_matched_cycles(BM, dim, index, plot_critical_values=False, use='mayavi'):
    vertices0, vertices1 = BM.get_matched_cycles(dim, index)

    if (len(vertices0) == 0 and len(vertices1) == 0):
        return

    if (len(BM.shape) == 2):
        cycle0 = np.zeros(shape=(BM.shape[0], BM.shape[1]))
        for vertex in vertices0[0:-1]:
            cycle0[vertex[0], vertex[1]] = 1
        if (plot_critical_values):
            if (len(vertices0) != 0):
                cycle0[vertices0[0][0], vertices0[0][1]] = 1.5
                cycle0[vertices0[-1][0], vertices0[-1][1]] = 0.5

        cycle1 = np.zeros(shape=(BM.shape[0], BM.shape[1]))
        for vertex in vertices1[0:-1]:
            cycle1[vertex[0], vertex[1]] = 1
        if (plot_critical_values):
            if (len(vertices1) != 0):
                cycle1[vertices1[0][0], vertices1[0][1]] = 1.5
                cycle1[vertices1[-1][0], vertices1[-1][1]] = 0.5
        
        fig = plt.figure(figsize=(15,5))
        fig.add_subplot(1, 2, 1)
        if (plot_critical_values):
            plt.imshow(cycle0, cmap='gray', vmin=0, vmax=1.5)
        else:
            plt.imshow(cycle0, cmap='gray', vmin=0, vmax=1)
        plt.axis('off')

        fig.add_subplot(1, 2, 2)
        if (plot_critical_values):
            plt.imshow(cycle1, cmap='gray', vmin=0, vmax=1.5)
        else:
            plt.imshow(cycle1, cmap='gray', vmin=0, vmax=1.5)
        plt.axis('off')

    if (len(BM.shape) ==3 ):
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

        if (use == 'mayavi'):
            cycle = np.concatenate([cycle0, cycle1], axis=0)
            fig = mlab.figure(size=(1000, 1000))
            contour = mlab.contour3d(cycle, colormap='bone', transparent=True)
            mlab.xlabel('X')
            mlab.ylabel('Y')
            mlab.zlabel('Z')
            mlab.show()

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
                surface_count=20,
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
                surface_count=20,
                caps= dict(x_show=True, y_show=True, z_show=True),
                ))
            fig = make_subplots(
            rows=1, cols=2,
            specs=[[{'type': 'volume'}, {'type': 'volume'}]])
            fig.add_trace(fig0.data[0], row=1, col=1)
            fig.add_trace(fig1.data[0], row=1, col=2)
            #fig.show()
            fig.write_html("plot.html")

    return