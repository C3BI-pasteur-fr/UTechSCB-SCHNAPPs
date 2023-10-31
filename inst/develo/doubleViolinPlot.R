library(plotly)

df <- read.csv("https://raw.githubusercontent.com/plotly/datasets/master/violin_data.csv")

pointposMale <- c(-0.3,-1.1,-0.6,-0.3)
pointposFemale <- c(0.45,0.55,1,0.4)
showLegend <- c(T,F,F,F)

fig <- plot_ly(type = 'violin')

i = 0
for (i in 1:length(unique(df$day))) {
  fig <- add_trace(
    fig,
    x = df$day[df$sex == 'Male' & df$day == unique(df$day)[i]],
    y = df$total_bill[df$sex == 'Male' & df$day == unique(df$day)[i]],
    hoveron = "points+kde",
    legendgroup = 'M',
    scalegroup = 'M',
    name = 'M',
    side = 'negative',
    box = list(
      visible = T
    ),
    points = 'all',
    pointpos = pointposMale[i],
    jitter = 0,
    scalemode = 'count',
    meanline = list(
      visible = T
    ),
    color = I("#8dd3c7"),
    marker = list(
      line = list(
        width = 2,
        color = "#8dd3c7"
      ),
      symbol = 'line-ns'
    ),
    showlegend = showLegend[i]
  ) 
  
  fig <- fig %>%
    add_trace(
      x = df$day[df$sex == 'Female' & df$day == unique(df$day)[i]],
      y = df$total_bill[df$sex == 'Female' & df$day == unique(df$day)[i]],
      hoveron = "points+kde",
      legendgroup = 'F',
      scalegroup = 'F',
      name = 'F',
      side = 'positive',
      box = list(
        visible = T
      ),
      points = 'all',
      pointpos = pointposFemale[i],
      jitter = 0,
      scalemode = 'count',
      meanline = list(
        visible = T
      ),
      color = I("#bebada"),
      marker = list(
        line = list(
          width = 2,
          color = "#bebada"
        ),
        symbol = 'line-ns'
      ),
      showlegend = showLegend[i]
    )
}

fig <- layout(
  fig,
  title = "Total bill distribution<br><i>scaled by number of bills per gender",
  yaxis = list(
    zeroline = F
  ),
  violingap = 0,
  violingroupgap = 0,
  violinmode = 'overlay',
  legend = list(
    tracegroupgap = 0
  )
)

fig

