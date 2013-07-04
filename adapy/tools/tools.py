import colorsys as _cs

def get_colors(color):
    for hue in range(color):
        hue = 1. * hue / color
        col = [int(x) for x in _cs.hsv_to_rgb(hue, 1.0, 230)]
        yield "#{0:02x}{1:02x}{2:02x}".format(*col)
    
