from vedo import Plotter, Points

def visualize_actuation(points):
    for i in range(len(points)):
        plt = Plotter(N=1,axes = 1)
        plt.show([Points(points[i])],f"time = {i}", at = 0)
        plt.close()