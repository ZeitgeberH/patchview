## base on https://stackoverflow.com/questions/48474699/marker-size-alpha-scaling-with-window-size-zoom-in-plot-scatter
import numpy as np
class FigureUpdater:
    def __init__(self):
        ##for storing information about Figures and Axes
        self.figs = {}
        ##for storing timers
        self.timer_dict = {}
        self.previous_scalingFactors = {}
    def add_ax(self, ax, features=[]):
        ax_dict = self.figs.setdefault(ax.figure,dict())
        ax_dict[ax] = {
            'xlim' : ax.get_xlim(),
            'ylim' : ax.get_ylim(),
            'figw' : ax.figure.get_figwidth(),
            'figh' : ax.figure.get_figheight(),
            'scale_s' : 1.0,
            'scale_a' : 1.0,
            'scale_linewidth' : 1.0,
            'features' : [features] if isinstance(features,str) else features,
        }
        self.previous_scalingFactors.update({ax:1.0})

    def update_axes(self):
        for fig,axes in self.figs.items():
            for ax, args in axes.items():
                fw = fig.get_figwidth()
                fh = fig.get_figheight()
                fac1 = min(fw/args['figw'], fh/args['figh'])
                xl = ax.get_xlim()
                yl = ax.get_ylim()
                fac2 = min(
                    abs(args['xlim'][1]-args['xlim'][0])/abs(xl[1]-xl[0]),
                    abs(args['ylim'][1]-args['ylim'][0])/abs(yl[1]-yl[0])
                )
                ##factor for marker size
                facS = (fac1*fac2)/args['scale_s']
                ##factor for alpha -- limited to values smaller 1.0
                facA = min(1.0,fac1*fac2)/args['scale_a']
                ## linewidth factor
                facL = (fac1*fac2)/args['scale_linewidth']
                ##updating the artists
                if facS != self.previous_scalingFactors[ax]:
                    for path in ax.collections:
                        if 'linewidth' in args['features']:
                            linewidth = path.get_linewidth()
                            if linewidth is not None:
                                path.set_linewidth(np.array(linewidth)*facL.tolist())
                    args['scale_s'] *= facS
                    args['scale_a'] *= facA
                    args['scale_linewidth'] *= facS
            self._redraw_later(fig)
            self.previous_scalingFactors[ax] = facS

    def _redraw_later(self, fig):
        timer = fig.canvas.new_timer(interval=10)
        timer.single_shot = True
        timer.add_callback(lambda : fig.canvas.draw_idle())
        timer.start()
        ##stopping previous timer
        if fig in self.timer_dict:
            self.timer_dict[fig].stop()
        ##storing a reference to prevent garbage collection
        self.timer_dict[fig] = timer