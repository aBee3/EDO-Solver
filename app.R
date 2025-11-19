##############################################
#     EDO SOLVER - FÁRMACO Y CAFEÍNA        #
##############################################
library(shiny)
library(shinythemes)
library(deSolve)
library(ggplot2)
library(gridExtra)

##############################################
# Solución Exacta - MODELO FÁRMACO
##############################################

Xg_exact <- function(t, a, X0, t0) {
  X0[1] * exp(-a * (t - t0))
}

Xb_exact <- function(t, a, k, X0, t0) {
  xg0 <- X0[1]
  xb0 <- X0[2]
  dt <- t - t0
  
  if (a == k) {
    return(xb0 * exp(-k * dt) + a * xg0 * dt * exp(-a * dt))
  } else {
    return(
      xb0 * exp(-k * dt) + 
        (a * xg0 / (k - a)) * (exp(-a * dt) - exp(-k * dt))
    )
  }
}

##############################################
# Solución Exacta - MODELO CAFEÍNA
##############################################

y_exact_cafeina <- function(t, k, y0, t0) {
  y0 * exp(-k * (t - t0))
}

##############################################
# Modelos EDO
##############################################

# Modelo Fármaco (sistema 2D)
m_farmaco <- function(a, k) {
  function(t, X) {
    Xg <- X[1]
    Xb <- X[2]
    dXg <- -a * Xg
    dXb <-  a * Xg - k * Xb
    return(c(dXg, dXb))
  }
}

# Modelo Cafeína (ecuación 1D)
m_cafeina <- function(k) {
  function(t, y) {
    return(-k * y)
  }
}

##############################################
# Métodos con deSolve
##############################################

# deSolve - Euler
desolve_euler <- function(f_ode, t0, X0, t_end, h) {
  times <- seq(t0, t_end, by = h)
  
  out <- ode(y = X0, times = times, func = function(t, y, parms) {
    list(f_ode(t, y))
  }, parms = NULL, method = "euler")  
  
  if (length(X0) == 1) {
    return(data.frame(step = 0:(nrow(out)-1), t = out[,1], y = out[,2]))
  } else {
    return(data.frame(step = 0:(nrow(out)-1), t = out[,1], Xg = out[,2], Xb = out[,3]))
  }
}

# deSolve - RK4
desolve_rk4 <- function(f_ode, t0, X0, t_end, h) {
  times <- seq(t0, t_end, by = h)
  
  out <- ode(y = X0, times = times, func = function(t, y, parms) {
    list(f_ode(t, y))
  }, parms = NULL, method = "rk4")  
  
  if (length(X0) == 1) {
    return(data.frame(step = 0:(nrow(out)-1), t = out[,1], y = out[,2]))
  } else {
    return(data.frame(step = 0:(nrow(out)-1), t = out[,1], Xg = out[,2], Xb = out[,3]))
  }
}

# deSolve - RK45 (ode45)
desolve_rk45 <- function(f_ode, t0, X0, t_end, h) {
  times <- seq(t0, t_end, by = h)
  
  out <- ode(y = X0, times = times, func = function(t, y, parms) {
    list(f_ode(t, y))
  }, parms = NULL, method = "ode45")  
  
  if (length(X0) == 1) {
    return(data.frame(step = 0:(nrow(out)-1), t = out[,1], y = out[,2]))
  } else {
    return(data.frame(step = 0:(nrow(out)-1), t = out[,1], Xg = out[,2], Xb = out[,3]))
  }
}

##############################################
# Métodos Numéricos
##############################################

# Euler 
euler_system <- function(f, t0, X0, h, n) {
  t <- numeric(n + 1)
  X <- matrix(0, n + 1, length(X0))
  
  t[1] <- t0
  X[1,] <- X0
  
  for (i in 1:n) {
    X[i+1,] <- X[i,] + h * f(t[i], X[i,])
    t[i+1] <- t[i] + h
  }
  
  if (length(X0) == 1) {
    return(data.frame(step=0:n, t=t, y=X[,1]))
  } else {
    return(data.frame(step=0:n, t=t, Xg=X[,1], Xb=X[,2]))
  }
}

# RK4 
rk4_system <- function(f, t0, X0, h, n) {
  t <- numeric(n + 1)
  X <- matrix(0, n + 1, length(X0))
  
  t[1] <- t0
  X[1,] <- X0
  
  for (i in 1:n) {
    k1 <- f(t[i], X[i,])
    k2 <- f(t[i] + h/2, X[i,] + h*k1/2)
    k3 <- f(t[i] + h/2, X[i,] + h*k2/2)
    k4 <- f(t[i] + h,   X[i,] + h*k3)
    
    X[i+1,] <- X[i,] + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
    t[i+1] <- t[i] + h
  }
  
  if (length(X0) == 1) {
    return(data.frame(step=0:n, t=t, y=X[,1]))
  } else {
    return(data.frame(step=0:n, t=t, Xg=X[,1], Xb=X[,2]))
  }
}

# RKDP45 
tolerance <- 1e-6
rk45_dp_system <- function(f, t0, X0, t_end, h_init, tol= tolerance) {
  
  a <- c(1/5, 3/10, 4/5, 8/9, 1, 1)
  b <- matrix(c(
    1/5,0,0,0,0,0,
    3/40,9/40,0,0,0,0,
    44/45,-56/15,32/9,0,0,0,
    19372/6561,-25360/2187,64448/6561,-212/729,0,0,
    9017/3168,-355/33,46732/5247,49/176,-5103/18656,0,
    35/384,0,500/1113,125/192,-2187/6784,11/84
  ), nrow=6, byrow=TRUE)
  
  c4 <- c(35/384,0,500/1113,125/192,-2187/6784,11/84,0)
  c5 <- c(5179/57600,0,7571/16695,393/640,-92097/339200,187/2100,1/40)
  
  t <- t0
  X <- X0
  h <- h_init
  
  if (length(X0) == 1) {
    sol <- data.frame(step=0, t=t, y=X[1], h=h)
  } else {
    sol <- data.frame(step=0, t=t, Xg=X[1], Xb=X[2], h=h)
  }
  
  while (t < t_end) {
    if (t + h > t_end) h <- t_end - t
    
    k1 <- f(t, X)
    k2 <- f(t + a[1]*h, X + h*(b[1,1]*k1))
    k3 <- f(t + a[2]*h, X + h*(b[2,1]*k1 + b[2,2]*k2))
    k4 <- f(t + a[3]*h, X + h*(b[3,1]*k1 + b[3,2]*k2 + b[3,3]*k3))
    k5 <- f(t + a[4]*h, X + h*(b[4,1]*k1 + b[4,2]*k2 + b[4,3]*k3 + b[4,4]*k4))
    k6 <- f(t + a[5]*h, X + h*(b[5,1]*k1 + b[5,2]*k2 + b[5,3]*k3 + b[5,4]*k4 + b[5,5]*k5))
    k7 <- f(t + a[6]*h, X + h*(b[6,1]*k1 + b[6,3]*k3 + b[6,4]*k4 + b[6,5]*k5 + b[6,6]*k6))
    
    X4 <- X + h*(c4[1]*k1 + c4[2]*k2 + c4[3]*k3 + c4[4]*k4 +
                   c4[5]*k5 + c4[6]*k6 + c4[7]*k7)
    
    X5 <- X + h*(c5[1]*k1 + c5[2]*k2 + c5[3]*k3 + c5[4]*k4 +
                   c5[5]*k5 + c5[6]*k6 + c5[7]*k7)
    
    err <- max(abs(X5 - X4))
    
    if (err < tol || h < 1e-12) {
      t <- t + h
      X <- X5
      if (length(X0) == 1) {
        sol <- rbind(sol, data.frame(step=nrow(sol), t=t, y=X[1], h=h))
      } else {
        sol <- rbind(sol, data.frame(step=nrow(sol), t=t, Xg=X[1], Xb=X[2], h=h))
      }
    }
    
    s <- 0.9 * (tol/(err + 1e-15))^(1/5)
    h <- h * max(0.1, min(5, s))
  }
  
  sol
}

##############################################
#               SHINY APP UI
##############################################

ui <- navbarPage(
  title = "EDO Solver: Euler Rk4 RK45 y desolve",
  theme = shinytheme("flatly"),
  
  tabPanel("Métodos",
           sidebarLayout(
             sidebarPanel(
               h3("Selección del Modelo"),
               selectInput("model", "Modelo:", 
                           choices = c("Fármaco", "Cafeína"), 
                           selected = "Fármaco"),
               
               hr(),
               h3("Parámetros del modelo"),
               
               conditionalPanel(
                 condition = "input.model == 'Fármaco'",
                 numericInput("a", "a (absorción)", value = 0.03, min = 0, max = 1, step = 0.001),
                 numericInput("k", "k (eliminación)", value = 0.02, min = 0, max = 1, step = 0.001),
                 numericInput("Xg0", "Xg(t0)", value = 400),
                 numericInput("Xb0", "Xb(t0)", value = 0)
               ),
               
               conditionalPanel(
                 condition = "input.model == 'Cafeína'",
                 numericInput("k_caf", "k (eliminación)", value = 0.15, min = 0, max = 1, step = 0.01),
                 numericInput("y0", "Concentración inicial (mg)", value = 70, min = 0)
               ),
               
               numericInput("t0", "Tiempo inicial", value = 0),
               numericInput("t_end", "Tiempo final", value = 350),
               numericInput("h", "Paso h", value = 20, min = 0.001, max = 50),
               
               actionButton("calcular", "CALCULAR", class="btn btn-primary btn-lg")
             ),
             
             mainPanel(
               h3("Comparación de Métodos Numéricos"),
               p("Presiona CALCULAR para ejecutar todos los métodos."),
               fluidRow(
                 column(width = 12,
                        plotOutput("plot_comparativo", height = "1200px", width = "100%")
                 )
               )
             )
           )
  ),
  
  tabPanel("Tabla de Pasos",
           h3("Tabla de pasos"),
           fluidRow(
             column(4,
                    h4("Euler"),
                    tableOutput("tabla_euler")
             ),
             column(4,
                    h4("RK4"),
                    tableOutput("tabla_rk4")
             ),
             column(4,
                    h4("RKDP45"),
                    tableOutput("tabla_rkdp45")
             )
           )
  ),
  
  tabPanel("Error Absoluto",   
           h3("Error absoluto"),
           fluidRow(
             column(12,
                    h4("Comparación del Error"),
                    plotOutput("plot_error_combined")
             )
           ),
           fluidRow(
             column(4,
                    h4("Error en Euler"),
                    plotOutput("plot_error_euler"),
                    h4("Tabla de pasos Euler"),
                    tableOutput("tabla_err_euler")
             ),
             column(4,
                    h4("Error en RK4"),
                    plotOutput("plot_error_rk4"),
                    h4("Tabla de pasos RK4"),
                    tableOutput("tabla_err_rk4")
             ),
             column(4,
                    h4("Error en RKDP45"),
                    plotOutput("plot_error_rk45"),
                    h4("Tabla de pasos Rk45"),
                    tableOutput("tabla_err_rk45")
             )
           )
  ),
  
  tabPanel("Learn More",
           fluidRow(
             h3("Visualización de la implementación del código"),
             column(5,
                    h4("Método de Euler"),
                    verbatimTextOutput("code_euler"),
                    h4("Método RK4"),
                    verbatimTextOutput("code_rk4")
             ),
             column(7,
                    h4("Método RKDP45"),
                    verbatimTextOutput("code_rkdp")
             )
           )
  ),  
  
  tabPanel("Conclusiones",
           fluidRow(
             h3("Conclusiones"),
             column(12,
                    h4("Comparativa de los métodos"),
                    p("El presente proyecto implementó y comparó tres métodos numéricos (Euler, Runge-Kutta 4 y Dormand-Prince RK45) para la solución del modelado de un farmaco y del metabolismo de la cafeina
                      y eliminación de cafeína. La comparación se realizó utilizando implementaciones propias de los algoritmos y las implementaciones estándar provistas por el paquete deSolve de R."),
                      br(),
                      p("Las simulaciones con paso fijo obtenidas mediante las funciones Euler y RK4 implementadas manualmente mostraron una concordancia casi perfecta con las soluciones generadas por las funciones de deSolve.
                      Esto confirma que estan desarrollados correctamente nuestros algoritmos."),
                      br(),
                      p("Método de Euler (Orden 1): Demostró ser el menos preciso, exhibiendo el mayor error absoluto, especialmente notable en los 
                      puntos críticos de las curvas (como el Cmax del modelo farmacológico). Un paso de tiempo grande (h=20) magnifica este error, limitando su utilidad para modelos donde la precisión es crucial."),
                      br(),
                      p("Método Runge-Kutta 4 (Orden 4): Proporcionó un excelente equilibrio entre complejidad y precisión. Su error absoluto es significativamente menor que el de Euler, ofreciendo soluciones virtualmente 
                      indistinguibles de la curva exacta en el rango de tiempo de interés para ambos modelos, incluso con un paso relativamente grande."),
                      br(),
                      p("Dormand-Prince RK45 (Adaptativo): Se reveló como el método más robusto y eficiente. Al ajustar el tamaño del paso (h) dinámicamente para mantener el error local dentro de la tolerancia (10^{-6}), minimiza el error global mientras optimiza el número de cálculos. Este
                      método es el recomendado para aplicaciones farmacocinéticas reales, donde la alta fidelidad en los resultados es indispensable."),
                      br(),
                      p("Para modelos dinámicos como el del metabolismo del fármaco, la elección del método numérico afecta directamente la validez de las predicciones. Si bien Euler es útil para la prueba de concepto, resulta inaceptable para el análisis cuantitativo. El método RK4 ofrece un buen compromiso 
                      de precisión con paso fijo. Sin embargo, el método RKDP45 (adaptativo) es el claro ganador en términos de precisión, robustez y eficiencia computacional, ya que garantiza la calidad del resultado sin requerir
                      que el usuario defina manualmente un paso óptimo."),
                  
             )
           )
  ),
  tabPanel("Individual",
           h3("Comparación Individual"),
           tabsetPanel(
             tabPanel("Euler",
                      fluidRow(
                        column(6,
                               h4("Solución vs Exacta"),
                               plotOutput("plot_euler_comparison", height = "600px")
                        ),
                        column(6,
                               h4("Error Absoluto"),
                               plotOutput("plot_euler_error", height = "600px")
                        )
                      ),
                      br(),
                      fluidRow(
                        column(12,
                               h4("Tabla de pasos y errores"),
                               tableOutput("table_euler_steps")
                        )
                      ),
                      br(),
                      fluidRow(
                        column(12,
                               h4("Código del método Euler"),
                               verbatimTextOutput("code_euler_copy")
                        )
                      )
             ),
             
             tabPanel("RK4",
                      fluidRow(
                        column(6,
                               h4("Solución vs Exacta"),
                               plotOutput("plot_rk4_comparison", height = "600px")
                        ),
                        column(6,
                               h4("Error Absoluto"),
                               plotOutput("plot_rk4_error", height = "600px")
                        )
                      ),
                      br(),
                      fluidRow(
                        column(12,
                               h4("Tabla de pasos y errores"),
                               tableOutput("table_rk4_steps")
                        )
                      ),
                      br(),
                      fluidRow(
                        column(12,
                               h4("Código del método Runge Kutta 4"),
                               verbatimTextOutput("code_rk4_copy")
                        )
                      )
             ),
             
             tabPanel("Dormand Prince Rk45",
                      fluidRow(
                        column(6,
                               h4("Solución vs Exacta"),
                               plotOutput("plot_rk45_comparison", height = "600px")
                        ),
                        column(6,
                               h4("Error Absoluto rk45"),
                               plotOutput("plot_rk45_error", height = "600px")
                        )
                      ),
                      br(),
                      fluidRow(
                        column(12,
                               h4("Código del método RK45"),
                               verbatimTextOutput("code_rk45_copy")
                        )
                      ),
                      br(),
                      fluidRow(
                        column(12,
                               h4("Tabla de pasos y errores"),
                               tableOutput("table_rk45_steps")
                        )
                      )
             )
           )
  )
)

##############################################
#              SERVER LOGIC
##############################################

server <- function(input, output, session) {
  
  values <- reactiveValues(
    df_euler = NULL,
    df_rk4 = NULL,
    df_rk45 = NULL,
    df_desolve_euler = NULL,
    df_desolve_rk4 = NULL,
    df_desolve_rk45 = NULL,
    err_euler = NULL,
    err_rk4 = NULL,
    err_rk45 = NULL,
    ready = FALSE
  )
  
  observeEvent(input$calcular, {
    
    req(input$h > 0)
    req(input$t_end > input$t0)
    
    isolate({
      
      t0 <- input$t0
      t_end <- input$t_end
      h <- input$h
      n <- as.integer((t_end - t0)/h)
      
      if (input$model == "Fármaco") {
        # MODELO FÁRMACO 
        a <- input$a
        k <- input$k
        X0 <- c(input$Xg0, input$Xb0)
        
        f <- m_farmaco(a, k)
        
        # Métodos propios
        values$df_euler <- euler_system(f, t0, X0, h, n)
        values$df_rk4   <- rk4_system(f, t0, X0, h, n)
        values$df_rk45  <- rk45_dp_system(f, t0, X0, t_end, h)
        
        # Métodos deSolve
        values$df_desolve_euler <- desolve_euler(f, t0, X0, t_end, h)
        values$df_desolve_rk4   <- desolve_rk4(f, t0, X0, t_end, h)
        values$df_desolve_rk45  <- desolve_rk45(f, t0, X0, t_end, h)
        
        # Calcular errores
        values$err_euler <- data.frame(
          step = seq_len(nrow(values$df_euler)),
          t = values$df_euler$t,
          xg_method = values$df_euler$Xg,
          xb_method = values$df_euler$Xb
        )
        values$err_euler$xg_exact <- Xg_exact(values$err_euler$t, a, X0, t0)
        values$err_euler$xb_exact <- Xb_exact(values$err_euler$t, a, k, X0, t0)
        values$err_euler$abs_err_xg <- abs(values$err_euler$xg_method - values$err_euler$xg_exact)
        values$err_euler$abs_err_xb <- abs(values$err_euler$xb_method - values$err_euler$xb_exact)
        
        values$err_rk4 <- data.frame(
          step = seq_len(nrow(values$df_rk4)),
          t = values$df_rk4$t,
          xg_method = values$df_rk4$Xg,
          xb_method = values$df_rk4$Xb
        )
        values$err_rk4$xg_exact <- Xg_exact(values$err_rk4$t, a, X0, t0)
        values$err_rk4$xb_exact <- Xb_exact(values$err_rk4$t, a, k, X0, t0)
        values$err_rk4$abs_err_xg <- abs(values$err_rk4$xg_method - values$err_rk4$xg_exact)
        values$err_rk4$abs_err_xb <- abs(values$err_rk4$xb_method - values$err_rk4$xb_exact)
        
        values$err_rk45 <- data.frame(
          step = values$df_rk45$step,
          t = values$df_rk45$t,
          xg_method = values$df_rk45$Xg,
          xb_method = values$df_rk45$Xb
        )
        values$err_rk45$xg_exact <- Xg_exact(values$err_rk45$t, a, X0, t0)
        values$err_rk45$xb_exact <- Xb_exact(values$err_rk45$t, a, k, X0, t0)
        values$err_rk45$abs_err_xg <- abs(values$err_rk45$xg_method - values$err_rk45$xg_exact)
        values$err_rk45$abs_err_xb <- abs(values$err_rk45$xb_method - values$err_rk45$xb_exact)
        
      } else {
        #  MODELO CAFEÍNA
        k_caf <- input$k_caf
        y0 <- input$y0
        X0 <- c(y0)
        
        f <- m_cafeina(k_caf)
        
        # Métodos propios
        values$df_euler <- euler_system(f, t0, X0, h, n)
        values$df_rk4   <- rk4_system(f, t0, X0, h, n)
        values$df_rk45  <- rk45_dp_system(f, t0, X0, t_end, h)
        
        # Métodos deSolve
        values$df_desolve_euler <- desolve_euler(f, t0, X0, t_end, h)
        values$df_desolve_rk4   <- desolve_rk4(f, t0, X0, t_end, h)
        values$df_desolve_rk45  <- desolve_rk45(f, t0, X0, t_end, h)
        
        # Calcular errores
        values$err_euler <- data.frame(
          step = seq_len(nrow(values$df_euler)),
          t = values$df_euler$t,
          y_method = values$df_euler$y
        )
        values$err_euler$y_exact <- y_exact_cafeina(values$err_euler$t, k_caf, y0, t0)
        values$err_euler$abs_err <- abs(values$err_euler$y_method - values$err_euler$y_exact)
        
        values$err_rk4 <- data.frame(
          step = seq_len(nrow(values$df_rk4)),
          t = values$df_rk4$t,
          y_method = values$df_rk4$y
        )
        values$err_rk4$y_exact <- y_exact_cafeina(values$err_rk4$t, k_caf, y0, t0)
        values$err_rk4$abs_err <- abs(values$err_rk4$y_method - values$err_rk4$y_exact)
        
        values$err_rk45 <- data.frame(
          step = values$df_rk45$step,
          t = values$df_rk45$t,
          y_method = values$df_rk45$y
        )
        values$err_rk45$y_exact <- y_exact_cafeina(values$err_rk45$t, k_caf, y0, t0)
        values$err_rk45$abs_err <- abs(values$err_rk45$y_method - values$err_rk45$y_exact)
      }
    })
    
    showNotification("Cálculo completado.", type = "message")
    values$ready <- TRUE
  })
  
  # Tablas de pasos
  output$tabla_euler <- renderTable({
    req(values$ready)   
    values$df_euler
  })
  
  output$tabla_rk4 <- renderTable({
    req(values$ready)
    values$df_rk4
  })
  
  output$tabla_rkdp45 <- renderTable({
    req(values$ready)
    values$df_rk45
  })
  
  # Plot comparativo principal
  output$plot_comparativo <- renderPlot({
    req(values$err_euler, values$err_rk4, values$err_rk45)
    req(values$df_desolve_euler, values$df_desolve_rk4, values$df_desolve_rk45)
    
    if (input$model == "Fármaco") {
      # GRÁFICA GRANDE COMBINADA 
      # Preparar datos para Xg
      df_combined_xg <- rbind(
        data.frame(t = values$err_euler$t, y = values$err_euler$xg_method, 
                   Metodo = "Euler", Tipo = "Implementación Propia", Compartimento = "Xg"),
        data.frame(t = values$df_desolve_euler$t, y = values$df_desolve_euler$Xg, 
                   Metodo = "Euler", Tipo = "deSolve", Compartimento = "Xg"),
        data.frame(t = values$err_rk4$t, y = values$err_rk4$xg_method, 
                   Metodo = "RK4", Tipo = "Implementación Propia", Compartimento = "Xg"),
        data.frame(t = values$df_desolve_rk4$t, y = values$df_desolve_rk4$Xg, 
                   Metodo = "RK4", Tipo = "deSolve", Compartimento = "Xg"),
        data.frame(t = values$err_rk45$t, y = values$err_rk45$xg_method, 
                   Metodo = "RK45", Tipo = "Implementación Propia", Compartimento = "Xg"),
        data.frame(t = values$df_desolve_rk45$t, y = values$df_desolve_rk45$Xg, 
                   Metodo = "RK45", Tipo = "deSolve", Compartimento = "Xg")
      )
      
      # Preparar datos para Xb
      df_combined_xb <- rbind(
        data.frame(t = values$err_euler$t, y = values$err_euler$xb_method, 
                   Metodo = "Euler", Tipo = "Implementación Propia", Compartimento = "Xb"),
        data.frame(t = values$df_desolve_euler$t, y = values$df_desolve_euler$Xb, 
                   Metodo = "Euler", Tipo = "deSolve", Compartimento = "Xb"),
        data.frame(t = values$err_rk4$t, y = values$err_rk4$xb_method, 
                   Metodo = "RK4", Tipo = "Implementación Propia", Compartimento = "Xb"),
        data.frame(t = values$df_desolve_rk4$t, y = values$df_desolve_rk4$Xb, 
                   Metodo = "RK4", Tipo = "deSolve", Compartimento = "Xb"),
        data.frame(t = values$err_rk45$t, y = values$err_rk45$xb_method, 
                   Metodo = "RK45", Tipo = "Implementación Propia", Compartimento = "Xb"),
        data.frame(t = values$df_desolve_rk45$t, y = values$df_desolve_rk45$Xb, 
                   Metodo = "RK45", Tipo = "deSolve", Compartimento = "Xb")
      )
      
      df_combined <- rbind(df_combined_xg, df_combined_xb)
      df_combined$Metodo_Completo <- paste0(df_combined$Metodo, " - ", df_combined$Compartimento)
      
      # Gráfica grande combinada
      p_combined <- ggplot(df_combined, aes(x = t, y = y, color = Metodo_Completo, 
                                            linetype = Tipo, size = Tipo)) +
        geom_line() +
        scale_color_manual(values = c(
          "Euler - Xg" = "#E41A1C", "Euler - Xb" = "#FF7F00",
          "RK4 - Xg" = "#377EB8", "RK4 - Xb" = "#4DAF4A",
          "RK45 - Xg" = "#984EA3", "RK45 - Xb" = "#A65628"
        )) +
        scale_linetype_manual(values = c("Implementación Propia" = "solid", "deSolve" = "dashed")) +
        scale_size_manual(values = c("Implementación Propia" = 1.2, "deSolve" = 0.8)) +
        labs(
          title = "Comparación Completa de Todos los Métodos Numéricos",
          subtitle = "Modelo Fármaco: Absorción Gastrointestinal y Sanguínea",
          x = "Tiempo (min)",
          y = "Concentración (mg)",
          color = "Método y Compartimento",
          linetype = "Implementación",
          size = "Implementación"
        ) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
          legend.position = "bottom",
          legend.box = "vertical",
          panel.grid.minor = element_line(color = "gray95"),
          panel.border = element_rect(color = "gray80", fill = NA, size = 0.5)
        ) +
        guides(
          color = guide_legend(nrow = 2, title.position = "top"),
          linetype = guide_legend(title.position = "top"),
          size = guide_legend(title.position = "top")
        )
      
      # 6 GRÁFICAS INDIVIDUALES
      # Euler Propio
      df_euler_propio <- data.frame(
        t = values$err_euler$t,
        Xg = values$err_euler$xg_method,
        Xb = values$err_euler$xb_method
      )
      df_euler_propio_long <- 
        data.frame(
          t = rep(df_euler_propio$t, times = ncol(df_euler_propio)-1),
          Compartimento = rep(names(df_euler_propio)[-1], each = nrow(df_euler_propio)),
          Concentracion = unlist(df_euler_propio[-1])
        )
      
      
      p1 <- ggplot(df_euler_propio_long, aes(x = t, y = Concentracion, color = Compartimento)) +
        geom_line(size = 1.2, linetype = "dashed") +
        scale_color_manual(values = c("Xg" = "#E41A1C", "Xb" = "#377EB8")) +
        labs(title = "Euler (Implementación Propia)", x = "Tiempo (min)", y = "Concentración (mg)") +
        theme_minimal(base_size = 11) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              legend.position = "top",
              panel.border = element_rect(color = "gray80", fill = NA, size = 0.5))
      
      # Euler deSolve
      df_euler_desolve <- data.frame(
        t = values$df_desolve_euler$t,
        Xg = values$df_desolve_euler$Xg,
        Xb = values$df_desolve_euler$Xb
      )
      df_euler_desolve_long <- data.frame(
        t = rep(df_euler_desolve$t, times = ncol(df_euler_desolve) - 1),
        Compartimento = rep(names(df_euler_desolve)[-1], each = nrow(df_euler_desolve)),
        Concentracion = unlist(df_euler_desolve[-1])
      )
      
      
      p2 <- ggplot(df_euler_desolve_long, aes(x = t, y = Concentracion, color = Compartimento)) +
        geom_line(size = 0.8, linetype = "dotted") +
        scale_color_manual(values = c("Xg" = "#E41A1C", "Xb" = "#377EB8")) +
        labs(title = "Euler (deSolve)", x = "Tiempo (min)", y = "Concentración (mg)") +
        theme_minimal(base_size = 11) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              legend.position = "top",
              panel.border = element_rect(color = "gray80", fill = NA, size = 0.5))
      
      # RK4 Propio
      df_rk4_propio <- data.frame(
        t = values$err_rk4$t,
        Xg = values$err_rk4$xg_method,
        Xb = values$err_rk4$xb_method
      )
      df_rk4_propio_long <- data.frame(
        t = rep(df_rk4_propio$t, times = ncol(df_rk4_propio) - 1),
        Compartimento = rep(names(df_rk4_propio)[-1], each = nrow(df_rk4_propio)),
        Concentracion = unlist(df_rk4_propio[-1])
      )
      
      
      p3 <- ggplot(df_rk4_propio_long, aes(x = t, y = Concentracion, color = Compartimento)) +
        geom_line(size = 1.2, linetype = "dashed") +
        scale_color_manual(values = c("Xg" = "#FF7F00", "Xb" = "#4DAF4A")) +
        labs(title = "RK4 (Implementación Propia)", x = "Tiempo (min)", y = "Concentración (mg)") +
        theme_minimal(base_size = 11) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              legend.position = "top",
              panel.border = element_rect(color = "gray80", fill = NA, size = 0.5))
      
      # RK4 deSolve
      df_rk4_desolve <- data.frame(
        t = values$df_desolve_rk4$t,
        Xg = values$df_desolve_rk4$Xg,
        Xb = values$df_desolve_rk4$Xb
      )
      df_rk4_desolve_long <- data.frame(
        t = rep(df_rk4_desolve$t, times = ncol(df_rk4_desolve) - 1),
        Compartimento = rep(names(df_rk4_desolve)[-1], each = nrow(df_rk4_desolve)),
        Concentracion = unlist(df_rk4_desolve[-1])
      )
      
      
      p4 <- ggplot(df_rk4_desolve_long, aes(x = t, y = Concentracion, color = Compartimento)) +
        geom_line(size = 0.8, linetype = "dotted") +
        scale_color_manual(values = c("Xg" = "#FF7F00", "Xb" = "#4DAF4A")) +
        labs(title = "RK4 (deSolve)", x = "Tiempo (min)", y = "Concentración (mg)") +
        theme_minimal(base_size = 11) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              legend.position = "top",
              panel.border = element_rect(color = "gray80", fill = NA, size = 0.5))
      
      # RK45 Propio
      df_rk45_propio <- data.frame(
        t = values$err_rk45$t,
        Xg = values$err_rk45$xg_method,
        Xb = values$err_rk45$xb_method
      )
      df_rk45_propio_long <- data.frame(
        t = rep(df_rk45_propio$t, times = ncol(df_rk45_propio) - 1),
        Compartimento = rep(names(df_rk45_propio)[-1], each = nrow(df_rk45_propio)),
        Concentracion = unlist(df_rk45_propio[-1])
      )
      
      p5 <- ggplot(df_rk45_propio_long, aes(x = t, y = Concentracion, color = Compartimento)) +
        geom_line(size = 1.2, linetype = "dashed") +
        scale_color_manual(values = c("Xg" = "#984EA3", "Xb" = "#A65628")) +
        labs(title = "RK45 (Implementación Propia)", x = "Tiempo (min)", y = "Concentración (mg)") +
        theme_minimal(base_size = 11) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              legend.position = "top",
              panel.border = element_rect(color = "gray80", fill = NA, size = 0.5))
      
      # RK45 deSolve
      df_rk45_desolve <- data.frame(
        t = values$df_desolve_rk45$t,
        Xg = values$df_desolve_rk45$Xg,
        Xb = values$df_desolve_rk45$Xb
      )
      df_rk45_desolve_long <- data.frame(
        t = rep(df_rk45_desolve$t, times = ncol(df_rk45_desolve) - 1),
        Compartimento = rep(names(df_rk45_desolve)[-1], each = nrow(df_rk45_desolve)),
        Concentracion = unlist(df_rk45_desolve[-1])
      )
      
      
      p6 <- ggplot(df_rk45_desolve_long, aes(x = t, y = Concentracion, color = Compartimento)) +
        geom_line(size = 0.8, linetype = "dotted") +
        scale_color_manual(values = c("Xg" = "#984EA3", "Xb" = "#A65628")) +
        labs(title = "RK45 (deSolve)", x = "Tiempo (min)", y = "Concentración (mg)") +
        theme_minimal(base_size = 11) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              legend.position = "top",
              panel.border = element_rect(color = "gray80", fill = NA, size = 0.5))
      
      # Combinar todas las gráficas
      grid.arrange(
        p_combined,
        arrangeGrob(p1, p2, p3, p4, p5, p6, ncol = 2),
        heights = c(1.5, 2)
      )
      
    } else {
      # MODELO CAFEÍNA 
      # Preparar datos combinados
      df_combined <- rbind(
        data.frame(t = values$err_euler$t, y = values$err_euler$y_method, 
                   Metodo = "Euler", Tipo = "Implementación Propia"),
        data.frame(t = values$df_desolve_euler$t, y = values$df_desolve_euler$y, 
                   Metodo = "Euler", Tipo = "deSolve"),
        data.frame(t = values$err_rk4$t, y = values$err_rk4$y_method, 
                   Metodo = "RK4", Tipo = "Implementación Propia"),
        data.frame(t = values$df_desolve_rk4$t, y = values$df_desolve_rk4$y, 
                   Metodo = "RK4", Tipo = "deSolve"),
        data.frame(t = values$err_rk45$t, y = values$err_rk45$y_method, 
                   Metodo = "RK45", Tipo = "Implementación Propia"),
        data.frame(t = values$df_desolve_rk45$t, y = values$df_desolve_rk45$y, 
                   Metodo = "RK45", Tipo = "deSolve")
      )
      
      # Gráfica grande combinada
      p_combined <- ggplot(df_combined, aes(x = t, y = y, color = Metodo, 
                                            linetype = Tipo, size = Tipo)) +
        geom_line() +
        scale_color_manual(values = c("Euler" = "#E41A1C", "RK4" = "#377EB8", "RK45" = "#4DAF4A")) +
        scale_linetype_manual(values = c("Implementación Propia" = "solid", "deSolve" = "dashed")) +
        scale_size_manual(values = c("Implementación Propia" = 1.2, "deSolve" = 0.8)) +
        labs(
          title = "Comparación Completa de Todos los Métodos Numéricos",
          subtitle = "Modelo de Eliminación de Cafeína",
          x = "Tiempo (horas)",
          y = "Concentración (mg)",
          color = "Método",
          linetype = "Implementación",
          size = "Implementación"
        ) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
          legend.position = "bottom",
          legend.box = "horizontal",
          panel.grid.minor = element_line(color = "gray95"),
          panel.border = element_rect(color = "gray80", fill = NA, size = 0.5)
        )
      
      # 6 GRÁFICAS INDIVIDUALES 
      # Euler Propio
      p1 <- ggplot(data.frame(t = values$err_euler$t, y = values$err_euler$y_method), 
                   aes(x = t, y = y)) +
        geom_line(color = "#E41A1C", size = 1.2, linetype = "dashed") +
        labs(title = "Euler (Implementación Propia)", x = "Tiempo (horas)", y = "Concentración (mg)") +
        theme_minimal(base_size = 11) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              panel.border = element_rect(color = "gray80", fill = NA, size = 0.5))
      
      # Euler deSolve
      p2 <- ggplot(data.frame(t = values$df_desolve_euler$t, y = values$df_desolve_euler$y), 
                   aes(x = t, y = y)) +
        geom_line(color = "#E41A1C", size = 0.8, linetype = "dotted") +
        labs(title = "Euler (deSolve)", x = "Tiempo (horas)", y = "Concentración (mg)") +
        theme_minimal(base_size = 11) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              panel.border = element_rect(color = "gray80", fill = NA, size = 0.5))
      
      # RK4 Propio
      p3 <- ggplot(data.frame(t = values$err_rk4$t, y = values$err_rk4$y_method), 
                   aes(x = t, y = y)) +
        geom_line(color = "#377EB8", size = 1.2, linetype = "dashed") +
        labs(title = "RK4 (Implementación Propia)", x = "Tiempo (horas)", y = "Concentración (mg)") +
        theme_minimal(base_size = 11) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              panel.border = element_rect(color = "gray80", fill = NA, size = 0.5))
      
      # RK4 deSolve
      p4 <- ggplot(data.frame(t = values$df_desolve_rk4$t, y = values$df_desolve_rk4$y), 
                   aes(x = t, y = y)) +
        geom_line(color = "#377EB8", size = 0.8, linetype = "dotted") +
        labs(title = "RK4 (deSolve)", x = "Tiempo (horas)", y = "Concentración (mg)") +
        theme_minimal(base_size = 11) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              panel.border = element_rect(color = "gray80", fill = NA, size = 0.5))
      
      # RK45 Propio
      p5 <- ggplot(data.frame(t = values$err_rk45$t, y = values$err_rk45$y_method), 
                   aes(x = t, y = y)) +
        geom_line(color = "#4DAF4A", size = 1.2, linetype = "dashed") +
        labs(title = "RK45 (Implementación Propia)", x = "Tiempo (horas)", y = "Concentración (mg)") +
        theme_minimal(base_size = 11) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              panel.border = element_rect(color = "gray80", fill = NA, size = 0.5))
      
      # RK45 deSolve
      p6 <- ggplot(data.frame(t = values$df_desolve_rk45$t, y = values$df_desolve_rk45$y), 
                   aes(x = t, y = y)) +
        geom_line(color = "#4DAF4A", size = 0.8, linetype = "dotted") +
        labs(title = "RK45 (deSolve)", x = "Tiempo (horas)", y = "Concentración (mg)") +
        theme_minimal(base_size = 11) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              panel.border = element_rect(color = "gray80", fill = NA, size = 0.5))
      
      # Combinar todas las gráficas
      grid.arrange(
        p_combined,
        arrangeGrob(p1, p2, p3, p4, p5, p6, ncol = 2),
        heights = c(1.5, 2)
      )
    }
  })
      
  # Plots de error individuales
  output$plot_error_euler <- renderPlot({
    req(values$err_euler)
    
    if (input$model == "Fármaco") {
      plot(values$err_euler$t, values$err_euler$abs_err_xg, type="l", col="orange", lwd=2,
           main="Error Método Euler", xlab="Tiempo (t)", ylab="Error absoluto")
      lines(values$err_euler$t, values$err_euler$abs_err_xb, col="blue", lwd=2)
      legend("topright", legend=c("Xg", "Xb"), col=c("orange", "blue"), lty=1)
    } else {
      plot(values$err_euler$t, values$err_euler$abs_err, type="l", col="orange", lwd=2,
           main="Error Método Euler", xlab="Tiempo (t)", ylab="Error absoluto")
    }
  })
  
  output$plot_error_rk4 <- renderPlot({
    req(values$err_rk4)
    
    if (input$model == "Fármaco") {
      plot(values$err_rk4$t, values$err_rk4$abs_err_xg, type="l", col="orange", lwd=2,
           main="Error Método RK4", xlab="Tiempo (t)", ylab="Error absoluto")
      lines(values$err_rk4$t, values$err_rk4$abs_err_xb, col="blue", lwd=2)
      legend("topright", legend=c("Xg", "Xb"), col=c("orange", "blue"), lty=1)
    } else {
      plot(values$err_rk4$t, values$err_rk4$abs_err, type="l", col="orange", lwd=2,
           main="Error Método RK4", xlab="Tiempo (t)", ylab="Error absoluto")
    }
  })
  
  output$plot_error_rk45 <- renderPlot({
    req(values$err_rk45)
    
    if (input$model == "Fármaco") {
      plot(values$err_rk45$t, values$err_rk45$abs_err_xg, type="l", col="orange", lwd=2,
           main="Error RKDP45", xlab="Tiempo (t)", ylab="Error absoluto")
      lines(values$err_rk45$t, values$err_rk45$abs_err_xb, col="blue", lwd=2)
      legend("topright", legend=c("Xg", "Xb"), col=c("orange", "blue"), lty=1)
    } else {
      plot(values$err_rk45$t, values$err_rk45$abs_err, type="l", col="orange", lwd=2,
           main="Error RKDP45", xlab="Tiempo (t)", ylab="Error absoluto")
    }
  })
  
  output$plot_error_combined <- renderPlot({
    req(values$err_euler, values$err_rk4, values$err_rk45)
    
    if (input$model == "Fármaco") {
      plot(values$err_euler$t, values$err_euler$abs_err_xg, type="l", col="red", lwd=2,
           xlab="Tiempo (t)", ylab="Error absoluto", main="Comparación del Error")
      lines(values$err_euler$t, values$err_euler$abs_err_xb, col=adjustcolor("red", 0.5), lwd=2)
      lines(values$err_rk4$t, values$err_rk4$abs_err_xg, col="blue", lwd=2)
      lines(values$err_rk4$t, values$err_rk4$abs_err_xb, col=adjustcolor("blue", 0.5), lwd=2)
      lines(values$err_rk45$t, values$err_rk45$abs_err_xg, col="green", lwd=2)
      lines(values$err_rk45$t, values$err_rk45$abs_err_xb, col=adjustcolor("green", 0.5), lwd=2)
      legend("topright",
             legend = c("Euler Xg", "Euler Xb", "RK4 Xg", "RK4 Xb", "RK45 Xg", "RK45 Xb"),
             col = c("red", adjustcolor("red", 0.5), "blue", adjustcolor("blue", 0.5),
                     "green", adjustcolor("green", 0.5)),
             lty = 1, lwd = 2, bty = "n")
    } else {
      plot(values$err_euler$t, values$err_euler$abs_err, type="l", col="red", lwd=2,
           xlab="Tiempo (t)", ylab="Error absoluto", main="Comparación del Error")
      lines(values$err_rk4$t, values$err_rk4$abs_err, col="blue", lwd=2)
      lines(values$err_rk45$t, values$err_rk45$abs_err, col="green", lwd=2)
      legend("topright",
             legend = c("Euler", "RK4", "RK45"),
             col = c("red", "blue", "green"),
             lty = 1, lwd = 2, bty = "n")
    }
  })
  
  # Tablas de error
  output$tabla_err_euler <- renderTable({
    req(values$err_euler)
    if (input$model == "Fármaco") {
      values$err_euler[, c("step","t","abs_err_xg","abs_err_xb")]
    } else {
      values$err_euler[, c("step","t","abs_err")]
    }
  })
  
  output$tabla_err_rk4 <- renderTable({
    req(values$err_rk4)
    if (input$model == "Fármaco") {
      values$err_rk4[, c("step","t","abs_err_xg","abs_err_xb")]
    } else {
      values$err_rk4[, c("step","t","abs_err")]
    }
  })
  
  output$tabla_err_rk45 <- renderTable({
    req(values$err_rk45)
    if (input$model == "Fármaco") {
      values$err_rk45[, c("step","t","abs_err_xg","abs_err_xb")]
    } else {
      values$err_rk45[, c("step","t","abs_err")]
    }
  })
  
  # Código de los métodos
  output$code_euler <- renderText({
    '
# MÉTODO DE EULER
# El método de Euler es el método numérico más básico para resolver EDOs
# Utiliza la aproximación: y(t + h) ≈ y(t) + h * f(t, y)
# donde f es la derivada dy/dt

euler_system <- function(f, t0, X0, h, n) {
  # Inicializar vectores para tiempo y solución
  t <- numeric(n + 1)  # Vector de tiempos
  X <- matrix(0, n + 1, length(X0))  # Matriz para las soluciones
  
  # Condiciones iniciales
  t[1] <- t0
  X[1, ] <- X0
  
  # Iteración del método de Euler
  for (i in 1:n) {
    # Fórmula de Euler: X_nuevo = X_actual + h * derivada
    X[i + 1, ] <- X[i, ] + h * f(t[i], X[i, ])
    t[i + 1] <- t[i] + h  # Avanzar el tiempo
  }
  
  # Retornar como data frame
  return(data.frame(step = 0:n, t = t, y = X))
}
    '
  })
  
  output$code_rk4 <- renderText({
    '
# MÉTODO RUNGE-KUTTA 4 (RK4)
# Es un método de orden 4, mucho más preciso que Euler
# Calcula 4 pendientes (k1, k2, k3, k4) en cada paso
# y las combina con pesos específicos para mayor precisión

rk4_system <- function(f, t0, X0, h, n) {
  # Inicializar vectores
  t <- numeric(n + 1)
  X <- matrix(0, n + 1, length(X0))
  
  # Condiciones iniciales
  t[1] <- t0
  X[1, ] <- X0
  
  # Iteración del método RK4
  for (i in 1:n) {
    # Calcular las 4 pendientes
    k1 <- f(t[i], X[i, ])                    # Pendiente al inicio
    k2 <- f(t[i] + h/2, X[i, ] + h/2 * k1)  # Pendiente en el punto medio (usando k1)
    k3 <- f(t[i] + h/2, X[i, ] + h/2 * k2)  # Pendiente en el punto medio (usando k2)
    k4 <- f(t[i] + h,   X[i, ] + h * k3)    # Pendiente al final (usando k3)
    
    # Combinar las pendientes con pesos: (k1 + 2k2 + 2k3 + k4) / 6
    X[i + 1, ] <- X[i, ] + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
    t[i + 1] <- t[i] + h
  }
  
  return(data.frame(step = 0:n, t = t, y = X))
}
    '
  })
  
  output$code_rkdp <- renderText({
    '
# MÉTODO DORMAND-PRINCE RK45 (RKDP45)
# Es un método adaptativo que ajusta automáticamente el tamaño del paso
# Calcula dos soluciones (orden 4 y orden 5) y compara el error
# Si el error es pequeño, acepta el paso; si no, lo reduce

rk45_dp_system <- function(f, t0, X0, t_end, h_init, tol=1e-6) {
  # Coeficientes de Dormand-Prince para calcular las pendientes
  a <- c(1/5, 3/10, 4/5, 8/9, 1, 1)
  
  # Matriz de coeficientes b para las combinaciones de pendientes
  b <- matrix(c(
    1/5,0,0,0,0,0,
    3/40,9/40,0,0,0,0,
    44/45,-56/15,32/9,0,0,0,
    19372/6561,-25360/2187,64448/6561,-212/729,0,0,
    9017/3168,-355/33,46732/5247,49/176,-5103/18656,0,
    35/384,0,500/1113,125/192,-2187/6784,11/84
  ), nrow=6, byrow=TRUE)
  
  # Coeficientes para solución de orden 4 y orden 5
  c4 <- c(35/384,0,500/1113,125/192,-2187/6784,11/84,0)
  c5 <- c(5179/57600,0,7571/16695,393/640,-92097/339200,187/2100,1/40)
  
  # Inicialización
  t <- t0
  X <- X0
  h <- h_init
  sol <- data.frame(step=0, t=t, y=X, h=h)
  
  # Integración adaptativa
  while (t < t_end) {
    if (t + h > t_end) h <- t_end - t  # Ajustar último paso
    
    # Calcular las 7 pendientes k1 a k7
    k1 <- f(t, X)
    k2 <- f(t + a[1]*h, X + h*(b[1,1]*k1))
    k3 <- f(t + a[2]*h, X + h*(b[2,1]*k1 + b[2,2]*k2))
    k4 <- f(t + a[3]*h, X + h*(b[3,1]*k1 + b[3,2]*k2 + b[3,3]*k3))
    k5 <- f(t + a[4]*h, X + h*(b[4,1]*k1 + b[4,2]*k2 + b[4,3]*k3 + b[4,4]*k4))
    k6 <- f(t + a[5]*h, X + h*(b[5,1]*k1 + b[5,2]*k2 + b[5,3]*k3 + b[5,4]*k4 + b[5,5]*k5))
    k7 <- f(t + a[6]*h, X + h*(b[6,1]*k1 + b[6,3]*k3 + b[6,4]*k4 + b[6,5]*k5 + b[6,6]*k6))
    
    # Calcular soluciones de orden 4 y 5
    X4 <- X + h*(c4[1]*k1 + c4[2]*k2 + c4[3]*k3 + c4[4]*k4 + c4[5]*k5 + c4[6]*k6 + c4[7]*k7)
    X5 <- X + h*(c5[1]*k1 + c5[2]*k2 + c5[3]*k3 + c5[4]*k4 + c5[5]*k5 + c5[6]*k6 + c5[7]*k7)
    
    # Estimar el error comparando ambas soluciones
    err <- max(abs(X5 - X4))
    
    # Si el error es aceptable, avanzar al siguiente paso
    if (err < tol || h < 1e-12) {
      t <- t + h
      X <- X5  # Usar la solución de orden 5 (más precisa)
      sol <- rbind(sol, data.frame(step=nrow(sol), t=t, y=X, h=h))
    }
    
    # Ajustar el tamaño del paso según el error
    # Si error pequeño → aumentar h; si error grande → reducir h
    s <- 0.9 * (tol/(err + 1e-15))^(1/5)
    h <- h * max(0.1, min(5, s))  # Limitar cambios entre 0.1x y 5x
  }
  
  return(sol)
}
  '
  })
  
  # Plots individuales - Euler
  output$plot_euler_comparison <- renderPlot({
    req(values$err_euler)
    df <- values$err_euler
    
    if (input$model == "Fármaco") {
      plot(df$t, df$xg_method, type="l", col="blue", lwd=3, lty=2,
           xlab="Tiempo (t)", ylab="Cantidad (mg)", main="Euler vs Exacta")
      lines(df$t, df$xg_exact, col="turquoise", lwd=2)
      lines(df$t, df$xb_method, col="magenta", lwd=3, lty=2)
      lines(df$t, df$xb_exact, col="pink", lwd=2)
      legend("topright", legend=c("Xg método", "Xg exacto", "Xb método", "Xb exacto"),
             col=c("blue", "turquoise", "magenta", "pink"), lty=c(2,1,2,1), lwd=c(3,2,3,2))
    } else {
      plot(df$t, df$y_method, type="l", col="blue", lwd=3, lty=2,
           xlab="Tiempo (t)", ylab="Concentración (mg)", main="Euler vs Exacta - Cafeína")
      lines(df$t, df$y_exact, col="turquoise", lwd=2)
      legend("topright", legend=c("Método", "Exacta"),
             col=c("blue", "turquoise"), lty=c(2,1), lwd=c(3,2))
    }
  })
  
  output$plot_euler_error <- renderPlot({
    req(values$err_euler)
    df <- values$err_euler
    
    if (input$model == "Fármaco") {
      plot(df$t, df$abs_err_xg, type="l", col="orange", lwd=2,
           xlab="Tiempo (t)", ylab="Error", main="Error absoluto (Euler)")
      lines(df$t, df$abs_err_xb, col="red", lwd=2)
      legend("topright", legend=c("Xg", "Xb"), col=c("orange", "red"), lty=1)
    } else {
      plot(df$t, df$abs_err, type="l", col="orange", lwd=2,
           xlab="Tiempo (t)", ylab="Error", main="Error absoluto (Euler)")
    }
  })
  
  output$table_euler_steps <- renderTable({
    req(values$df_euler)
    head(values$df_euler, 80)
  })
  
  output$code_euler_copy <- renderText({
    '
    euler_system <- function(f, t0, X0, h, n) {
      t <- numeric(n + 1)
      X <- matrix(0, n + 1, length(X0))
      t[1] <- t0
      X[1, ] <- X0
      
      for (i in 1:n) {
        X[i + 1, ] <- X[i, ] + h * f(t[i], X[i, ])
        t[i + 1] <- t[i] + h
      }
      
      return(data.frame(step = 0:n, t = t, y = X))
    }
    '
  })
  
  # Plots individuales - RK4
  output$plot_rk4_comparison <- renderPlot({
    req(values$err_rk4)
    df <- values$err_rk4
    
    if (input$model == "Fármaco") {
      plot(df$t, df$xg_method, type="l", col="blue", lwd=4, lty=2,
           xlab="Tiempo (t)", ylab="Cantidad (mg)", main="RK4 vs Exacta")
      lines(df$t, df$xg_exact, col="turquoise", lwd=2)
      lines(df$t, df$xb_method, col="magenta", lwd=4, lty=2)
      lines(df$t, df$xb_exact, col="pink", lwd=2)
      legend("topright", legend=c("Xg método", "Xg exacto", "Xb método", "Xb exacto"),
             col=c("blue", "turquoise", "magenta", "pink"), lty=c(2,1,2,1), lwd=c(4,2,4,2))
    } else {
      plot(df$t, df$y_method, type="l", col="blue", lwd=4, lty=2,
           xlab="Tiempo (t)", ylab="Concentración (mg)", main="RK4 vs Exacta - Cafeína")
      lines(df$t, df$y_exact, col="turquoise", lwd=2)
      legend("topright", legend=c("Método", "Exacta"),
             col=c("blue", "turquoise"), lty=c(2,1), lwd=c(4,2))
    }
  })
  
  output$plot_rk4_error <- renderPlot({
    req(values$err_rk4)
    df <- values$err_rk4
    
    if (input$model == "Fármaco") {
      plot(df$t, df$abs_err_xg, type="l", col="orange", lwd=2,
           xlab="Tiempo (t)", ylab="Error", main="Error absoluto (RK4)")
      lines(df$t, df$abs_err_xb, col="red", lwd=2)
      legend("topright", legend=c("Xg", "Xb"), col=c("orange", "red"), lty=1)
    } else {
      plot(df$t, df$abs_err, type="l", col="orange", lwd=2,
           xlab="Tiempo (t)", ylab="Error", main="Error absoluto (RK4)")
    }
  })
  
  output$table_rk4_steps <- renderTable({
    req(values$df_rk4)
    head(values$df_rk4, 80)
  })
  
  output$code_rk4_copy <- renderText({
    '
    rk4_system <- function(f, t0, X0, h, n) {
      t <- numeric(n + 1)
      X <- matrix(0, n + 1, length(X0))
      t[1] <- t0
      X[1, ] <- X0
      
      for (i in 1:n) {
        k1 <- f(t[i], X[i, ])
        k2 <- f(t[i] + h/2, X[i, ] + h/2 * k1)
        k3 <- f(t[i] + h/2, X[i, ] + h/2 * k2)
        k4 <- f(t[i] + h,   X[i, ] + h * k3)
        
        X[i + 1, ] <- X[i, ] + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
        t[i + 1] <- t[i] + h
      }
      
      return(data.frame(step = 0:n, t = t, y = X))
    }
    '
  })
  
  # Plots individuales - RK45
  output$plot_rk45_comparison <- renderPlot({
    req(values$err_rk45)
    df <- values$err_rk45
    
    if (input$model == "Fármaco") {
      plot(df$t, df$xg_method, type="l", col="blue", lwd=4, lty=2,
           xlab="Tiempo (t)", ylab="Cantidad (mg)", main="RK45 vs Exacta")
      lines(df$t, df$xg_exact, col="turquoise", lwd=2)
      lines(df$t, df$xb_method, col="magenta", lwd=4, lty=2)
      lines(df$t, df$xb_exact, col="pink", lwd=2)
      legend("topright", legend=c("Xg método", "Xg exacto", "Xb método", "Xb exacto"),
             col=c("blue", "turquoise", "magenta", "pink"), lty=c(2,1,2,1), lwd=c(4,2,4,2))
    } else {
      plot(df$t, df$y_method, type="l", col="blue", lwd=4, lty=2,
           xlab="Tiempo (t)", ylab="Concentración (mg)", main="RK45 vs Exacta - Cafeína")
      lines(df$t, df$y_exact, col="turquoise", lwd=2)
      legend("topright", legend=c("Método", "Exacta"),
             col=c("blue", "turquoise"), lty=c(2,1), lwd=c(4,2))
    }
  })
  
  output$plot_rk45_error <- renderPlot({
    req(values$err_rk45)
    df <- values$err_rk45
    
    if (input$model == "Fármaco") {
      plot(df$t, df$abs_err_xg, type="l", col="orange", lwd=2,
           xlab="Tiempo (t)", ylab="Error", main="Error absoluto (RK45)")
      lines(df$t, df$abs_err_xb, col="red", lwd=2)
      legend("topright", legend=c("Xg", "Xb"), col=c("orange", "red"), lty=1)
    } else {
      plot(df$t, df$abs_err, type="l", col="orange", lwd=2,
           xlab="Tiempo (t)", ylab="Error", main="Error absoluto (RK45)")
    }
  })
  
  output$table_rk45_steps <- renderTable({
    req(values$df_rk45)
    head(values$df_rk45, 75)
  })
  
  output$code_rk45_copy <- renderText({
    '
  rk45_dp_system <- function(f, t0, X0, t_end, h_init, tol=1e-6) {
    a <- c(1/5, 3/10, 4/5, 8/9, 1, 1)
    b <- matrix(c(
      1/5,0,0,0,0,0,
      3/40,9/40,0,0,0,0,
      44/45,-56/15,32/9,0,0,0,
      19372/6561,-25360/2187,64448/6561,-212/729,0,0,
      9017/3168,-355/33,46732/5247,49/176,-5103/18656,0,
      35/384,0,500/1113,125/192,-2187/6784,11/84
    ), nrow=6, byrow=TRUE)
    
    c4 <- c(35/384,0,500/1113,125/192,-2187/6784,11/84,0)
    c5 <- c(5179/57600,0,7571/16695,393/640,-92097/339200,187/2100,1/40)
    
    t <- t0; X <- X0; h <- h_init
    sol <- data.frame(step=0, t=t, y=X, h=h)
    
    while (t < t_end) {
      if (t + h > t_end) h <- t_end - t
      
      k1 <- f(t, X)
      k2 <- f(t + a[1]*h, X + h*(b[1,1]*k1))
      k3 <- f(t + a[2]*h, X + h*(b[2,1]*k1 + b[2,2]*k2))
      k4 <- f(t + a[3]*h, X + h*(b[3,1]*k1 + b[3,2]*k2 + b[3,3]*k3))
      k5 <- f(t + a[4]*h, X + h*(b[4,1]*k1 + b[4,2]*k2 + b[4,3]*k3 + b[4,4]*k4))
      k6 <- f(t + a[5]*h, X + h*(b[5,1]*k1 + b[5,2]*k2 + b[5,3]*k3 + b[5,4]*k4 + b[5,5]*k5))
      k7 <- f(t + a[6]*h, X + h*(b[6,1]*k1 + b[6,3]*k3 + b[6,4]*k4 + b[6,5]*k5 + b[6,6]*k6))
      
      X4 <- X + h*(c4[1]*k1 + c4[2]*k2 + c4[3]*k3 + c4[4]*k4 + c4[5]*k5 + c4[6]*k6 + c4[7]*k7)
      X5 <- X + h*(c5[1]*k1 + c5[2]*k2 + c5[3]*k3 + c5[4]*k4 + c5[5]*k5 + c5[6]*k6 + c5[7]*k7)
      
      err <- max(abs(X5 - X4))
      
      if (err < tol || h < 1e-12) {
        t <- t + h; X <- X5
        sol <- rbind(sol, data.frame(step=nrow(sol), t=t, y=X, h=h))
      }
      
      s <- 0.9 * (tol/(err + 1e-15))^(1/5)
      h <- h * max(0.1, min(5, s))
    }
    return(sol)
  }
  '
  })
}

##############################################
# RUN APP
##############################################

shinyApp(ui, server)