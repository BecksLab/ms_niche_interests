library(dplyr)
library(tidyr)
library(purrr)
library(plotly)
library(colorspace)

pc_vars <- paste0("PC", 1:3)

##------------------------------------------------------------
## Ellipsoid surface
##------------------------------------------------------------

make_ellipsoid <- function(covmat, centre, n = 40, scale = 2){
  
  eig <- eigen(covmat)
  
  radii <- scale * sqrt(pmax(eig$values, 0))
  
  u <- seq(0, 2*pi, length.out = n)
  v <- seq(0, pi, length.out = n)
  
  x <- outer(cos(u), sin(v))
  y <- outer(sin(u), sin(v))
  z <- outer(rep(1, n), cos(v))
  
  pts <- cbind(c(x), c(y), c(z))
  
  pts <- pts %*% diag(radii)
  pts <- pts %*% t(eig$vectors)
  pts <- sweep(pts, 2, centre, "+")
  
  list(
    x = matrix(pts[,1], n, n),
    y = matrix(pts[,2], n, n),
    z = matrix(pts[,3], n, n)
  )
  
}

##------------------------------------------------------------
## Overall ellipsoids
##------------------------------------------------------------

state_summary <-
  topology_space %>%
  squad_up(state) %>%
  no_cap(
    centre = list(as.numeric(colMeans(across(all_of(pc_vars))))),
    covmat = list(cov(across(all_of(pc_vars)))),
    .groups = "drop"
  ) %>%
  glow_up(
    colour = c(
      Pre  = shark_silver,
      Post = darken(shark_silver, .30)
    )[state]
  )

##------------------------------------------------------------
## Model-specific ellipsoids
##------------------------------------------------------------

model_state_summary <-
  topology_space %>%
  squad_up(model, state) %>%
  no_cap(
    centre = list(as.numeric(colMeans(across(all_of(pc_vars))))),
    covmat = list(cov(across(all_of(pc_vars)))),
    .groups = "drop"
  ) %>%
  glow_up(
    colour = if_else(
      state == "Post",
      darken(model_colours[model], .30),
      model_colours[model]
    )
  )

##------------------------------------------------------------
## Model centroids
##------------------------------------------------------------

model_summary <-
  topology_space %>%
  squad_up(model, state) %>%
  no_cap(
    centre = list(as.numeric(colMeans(across(all_of(pc_vars))))),
    .groups = "drop"
  ) %>%
  glow_up(
    colour = if_else(
      state == "Post",
      darken(model_colours[model], .30),
      model_colours[model]
    )
  )

arrows <-
  model_summary %>%
  pivot_wider(
    names_from = state,
    values_from = c(centre, colour)
  )

##------------------------------------------------------------
## Generic plotting function
##------------------------------------------------------------

plot_scene <- function(ellipsoids,
                       arrows,
                       title = NULL){
  
  fig <- plot_ly()
  
  ##------------------------
  ## Ellipsoid surfaces
  ##------------------------
  
  for(i in seq_len(nrow(ellipsoids))){
    
    ell <- make_ellipsoid(
      ellipsoids$covmat[[i]],
      ellipsoids$centre[[i]]
    )
    
    fig <-
      fig %>%
      add_surface(
        x = ell$x,
        y = ell$y,
        z = ell$z,
        opacity = 0.18,
        showscale = FALSE,
        surfacecolor = matrix(0, nrow(ell$x), ncol(ell$x)),
        colorscale = list(
          c(0, ellipsoids$colour[i]),
          c(1, ellipsoids$colour[i])
        ),
        hoverinfo = "skip",
        showlegend = FALSE
      )
    
  }
  
  ##------------------------
  ## Model arrows
  ##------------------------
  
  for(i in seq_len(nrow(arrows))){
    
    pre  <- arrows$centre_Pre[[i]]
    post <- arrows$centre_Post[[i]]
    
    fig <-
      fig %>%
      add_trace(
        x = c(pre[1], post[1]),
        y = c(pre[2], post[2]),
        z = c(pre[3], post[3]),
        type = "scatter3d",
        mode = "lines+markers",
        line = list(
          color = arrows$colour_Pre[[i]],
          width = 8
        ),
        marker = list(
          size = 5,
          color = c(
            arrows$colour_Pre[[i]],
            arrows$colour_Post[[i]]
          )
        ),
        name = arrows$model[i]
      )
    
  }
  
  fig %>%
    layout(
      title = title,
      scene = list(
        xaxis = list(title = "PC1"),
        yaxis = list(title = "PC2"),
        zaxis = list(title = "PC3"),
        aspectmode = "data"
      )
    )
  
}

##------------------------------------------------------------
## Overall plot
##------------------------------------------------------------

overall_fig <-
  plot_scene(
    ellipsoids = state_summary,
    arrows = arrows,
    title = "Overall"
  )

overall_fig

##------------------------------------------------------------
## Individual model plots
##------------------------------------------------------------

models <- unique(topology_space$model)

model_plots <-
  map(models, function(m){
    
    plot_scene(
      ellipsoids = filter(model_state_summary, model == m),
      arrows     = filter(arrows, model == m),
      title      = m
    )
    
  })

names(model_plots) <- models

## Example:
model_plots[[4]]

##------------------------------------------------------------
## Multi-panel figure
##------------------------------------------------------------

## Build each widget into traces
plots <- lapply(model_plots, plotly_build)

models <- names(model_plots)

ncols <- 2
nrows <- ceiling(length(plots) / ncols)

fig <- plot_ly()

## Layout information
layout_list <- list()
annotations <- list()

for(i in seq_along(plots)) {
  
  scene <- if(i == 1) "scene" else paste0("scene", i)
  
  ##--------------------------------------------------
  ## Add every trace from this plot into the scene
  ##--------------------------------------------------
  
  for(tr in plots[[i]]$x$data) {
    
    tr$scene <- scene
    
    fig$x$data <- append(fig$x$data, list(tr))
    
  }
  
  ##--------------------------------------------------
  ## Position of this scene
  ##--------------------------------------------------
  
  row <- ceiling(i / ncols)
  col <- ((i - 1) %% ncols) + 1
  
  x0 <- (col - 1) / ncols
  x1 <- col / ncols
  
  ## Top row first
  y1 <- 1 - (row - 1) / nrows
  y0 <- 1 - row / nrows
  
  layout_list[[scene]] <- list(
    
    domain = list(
      x = c(x0, x1),
      y = c(y0, y1)
    ),
    
    xaxis = list(title = "PC1"),
    
    yaxis = list(title = "PC2"),
    
    zaxis = list(title = "PC3"),
    
    aspectmode = "data"
    
  )
  
  annotations[[i]] <- list(
    
    x = mean(c(x0, x1)),
    y = y1 - 0.01,
    
    xref = "paper",
    yref = "paper",
    
    text = paste0("<b>", models[i], "</b>"),
    
    showarrow = FALSE,
    
    font = list(size = 16)
    
  )
  
}

##--------------------------------------------------
## Inject layouts
##--------------------------------------------------

fig$x$layout <- c(
  fig$x$layout,
  layout_list,
  list(
    annotations = annotations,
    showlegend = FALSE
  )
)

fig
