
_INCHES_PER_PT = 1.0/72.27                               # convert pt to inches
_INCHES_PER_MM = 1/25.4                                  # convert mm to inches
_GOLDEN_MEAN = (5**0.5-1.0)/2.0   


class AandA(): pass
    
# ASTRONOMY & ASTROPHYSICS definitions for figure sizes
# 255.76535 pt column width equals 88 mm (from aa.doc.pdf section 3.4)
# one inch in mm = 25.4 mm/inch
#one_col_fig_width_pt = 249.448819
AandA.ONE_COL_FIG_WIDTH_MM = 88
AandA.TWO_COL_FIG_WIDTH_MM = 180
AandA.SIDE_CAPTION_FIG_WIDTH_MM = 120
AandA.ONE_COL_FIG_WIDTH = AandA.ONE_COL_FIG_WIDTH_MM * _INCHES_PER_MM  # width in inches
AandA.ONE_COL_FIG_HEIGHT = AandA.ONE_COL_FIG_WIDTH * _GOLDEN_MEAN      # height in inches
AandA.TWO_COL_FIG_WIDTH = AandA.TWO_COL_FIG_WIDTH_MM * _INCHES_PER_MM
AandA.SIDE_CAPTION_FIG_WIDTH = AandA.SIDE_CAPTION_FIG_WIDTH_MM * _INCHES_PER_MM
AandA.FIG_SIZE = [AandA.ONE_COL_FIG_WIDTH, AandA.ONE_COL_FIG_HEIGHT]

