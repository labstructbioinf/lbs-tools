from typing import Union, List, Tuple

import numpy as np
import matplotlib.pyplot as plt


class DivColorScaling:
    '''
    base class for coloring residue sequences according to their importance
    '''
    tdseq = '<td style="text-align:left; font-weight:bold; letter-spacing:1px;">'
    #table_style = '<table style="border-spacing: 10px, border-collapse: collapse;">'
    table_style = '<table>'
    cmap = plt.cm.Greens
    #cmap = sns.color_palette("coolwarm")
    #cmap = sns.diverging_palette(240, 140, s=80, l=55, n=5, center='dark', as_cmap=True)
    @staticmethod
    def scale_color(color):
        '''
        scale given number's to string ing 0-255
        '''
        if isinstance(color, (float, np.floating)):
            color = int(color*255)
        else:
            color=int(color)
        return str(color)
    
    def pick_color(self, imp):
        '''
        maps residue importance/score to RGB color sting
        Args:
            imp (float, int)
        Returns:
            rbg_color: (str) color string
        '''
        imp = self.cmap(imp)
        if len(imp)==4: #removes alpha channel from rgba tuple
            imp = imp[:3]
        imp = [self.scale_color(i) for i in imp]
        return 'rgb(' + ', '.join(imp) + ')'

    def html_colored_letter(self, letter: str, importance: Union[float, int]) -> str:
        '''
        add html styling for given aa letter
        Args:
            letter: (str)
            importance: (float, int) single number in range 0-1 for floats and 0-255 to ints
        Returns:
            html_string: (str)
        '''
        letter_color = self.pick_color(importance)
        html_string = '<span style= "color:' + letter_color + '">' + letter + '</span>'
        return html_string

    def scale_embeddings(self, emb, scale_factor=0.95):
        '''
        sign conserved embedding scaling
        params:
            scale_factor (float) 0-1.0 lower values increases cmap diversity
        '''
        if emb.shape[0] > 4:
            emb = emb.reshape(4, -1)
        sign = 1/(1 + np.exp(-emb.sum(0))).max()
        emb_abs_max = np.abs(emb).max()*sign*scale_factor
        x = 0.5 + 0.5*np.clip(emb/emb_abs_max, a_min=-1, a_max=1)
        return x

    def tdi(self, string):
        string_len = len(string)
        diff_len = 8 - string_len
        string_justed = string
        if diff_len > 0:
            string_justed += '&nbsp;'*diff_len
        return '<td>' + string_justed + '</td>'
