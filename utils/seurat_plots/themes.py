import plotnine as p9

umap_theme = p9.theme(
    panel_background=p9.element_rect(fill="white"),
    axis_line =p9.element_blank(),
    axis_text = p9.element_blank(),
    axis_ticks = p9.element_blank(),
    axis_title_x=p9.element_blank(),
    axis_title_y=p9.element_blank(),

)