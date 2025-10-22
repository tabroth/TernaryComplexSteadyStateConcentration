# Ternary complex steady state concentration ShinyApp
library(plotly)
library(shinyjs)

ui <- fluidPage(
  useShinyjs(),
  titlePanel("Ternary complex simulations"),
  withMathJax(),
  sidebarLayout(
    sidebarPanel(width = 3,
      tabsetPanel(
       tabPanel("Data", br(),
                 textAreaInput("Ligtot", "Titration series of ligand L (in nM):", value = "1\n5\n10\n25\n50\n100\n250\n500\n750\n1000\n2500\n5000\n7500\n10000\n25000\n50000\n75000\n100000"), 
                 textAreaInput("Tartot", "Total target concentrations T (in nM):", value = "50"), 
                 textAreaInput("E3tot", "Total chaperone concentrations C (in nM):", value = "50"), 
        ),
        tabPanel("Parameters", br(),
                 numericInput("KP1", "Equilibrium binding constant of ligand and target protein \\(K_\\text{T}\\) (in nM):", min = 0, 10000, step=1e0),
                 radioButtons("KT_choice", "Select target-ligand \\(K_\\text{D}\\):", selected = "KT1", 
                              choices = c(
                                "\\(K_\\text{T,1}\\)" = "KT1",
                                "\\(K_\\text{T,2}\\)=\\(K_\\text{T,1}/\\alpha\\)" = "KT2"),
                              inline = TRUE),
                 numericInput("KE1", "Equilibrium binding constant of ligand and chaperone \\(K_\\text{C,1}\\) (in nM):", min = 0, 100, step=1e1),
                 numericInput("alpha", "Cooperativity \\(\\alpha\\):", 1, min = 0, max = 300, step=1e0),
                 radioButtons("Native_interaction", "Native interaction between target and chaperone:", selected = "none", 
                              choices = c(
                                "None" = "none",
                                "Native" = "native")),    
                 numericInput("KPE", "Equilibrium dissociation constant of target and chaperone \\(K_\\text{TC} \\) (in nM):", min = 0, 1000, step=1e1),
                 radioButtons("Cooperativity_choice", "Select cooperativity:", selected = "alpha", 
                              choices = c(
                                "\\(\\alpha\\)" = "alpha",
                                "\\(\\beta\\)" = "beta",
                                "\\(\\gamma\\)" = "gamma"),
                              inline = TRUE),    
                 numericInput("Cooperativity", "Cooperativity:", min = 0, 1, step=1e1),
                 radioButtons("Tracer_binding", "Additional binding to chaperone:", selected = "none", 
                              choices = c(
                                "None" = "none",
                                "Competitive" = "competitive",
                                "Conformational equilibrium" = "conf_eq")),    
                 numericInput("KC", "Equilibrium dissociation constant of the competitor Tr for the chaperone \\(K_\\text{Tr} \\) (in nM):", min = 0, 65, step=1e1),
                 numericInput("Ctot", "Total concentration of tracer \\([\\text{Tr}_\\text{tot}] \\) (in nM):", min = 0, 25, step=1e0),
                 numericInput("Kconf", "Conformational equilibrium constant of chaperone and tracer \\(K_\\text{c} = [\\text{C}_\\text{c}]/[\\text{C}] \\):", min = 0, 100, step=1e1)                 ),
       tabPanel("Display", br(),
                radioButtons("Ligand_display", "Plot scale for x-axis", selected = "Total_ligand", 
                               choices = c(
                               "Free ligand" = "Free_ligand",
                               "Total ligand" = "Total_ligand",
                               "Total target" = "Total_target",
                               "Total chaperone" = "Total_chap")),    
       )
      )
    ),
    mainPanel(
      navbarPage(br(),
        tabPanel("Model scheme", align = "center", 
          img(src = "Scheme_Ternary.png", alt = "Scheme", width = 0.35*2382 , height = 0.35*1320)
        ),
        tabPanel("Complexes",
          tabsetPanel(type="tabs",
            tabPanel("Combined", align = "center", plotlyOutput("plot_cmplx")),
            tabPanel("Ternary", align = "center", plotlyOutput("plot_PLE")),
            tabPanel("Chaperone–ligand", align = "center", plotlyOutput("plot_EL")),
            tabPanel("Target–ligand", align = "center", plotlyOutput("plot_PL")),
            tabPanel("Target–chaperone", align = "center", plotlyOutput("plot_PE")),
            tabPanel("Ternary & target–ligand", align = "center", plotlyOutput("plot_PLE_PL")),
            tabPanel("Ternary & chaperone–ligand", align = "center", plotlyOutput("plot_PLE_EL")),
            tabPanel("Ternary & target–chaperone", align = "center", plotlyOutput("plot_PLE_PE"))
          )
        ),
        tabPanel("Target",
          tabsetPanel(type="tabs",
            tabPanel("Total", align = "center", plotlyOutput("plot_tar_tot")),
            tabPanel("Free", align = "center", plotlyOutput("plot_tar_free")),
            tabPanel("% bound", align = "center", plotlyOutput("plot_TO")),
            tabPanel("% target and chaperone bound", align = "center", plotlyOutput("plot_Tbnd_Cbnd")),
            tabPanel("% ternary complex", align = "center", plotlyOutput("plot_PLE_P"))
          )
        ),
        tabPanel("Ligand",
          tabsetPanel(type="tabs",
              tabPanel("Free-to-total conversion", align = "center", plotlyOutput("plot_lig_tot")),
              tabPanel("% bound", align = "center", plotlyOutput("plot_lig_bound"))
          )
        ),
        tabPanel("Chaperone",
                 tabsetPanel(type="tabs",
                             tabPanel("Total", align = "center", plotlyOutput("plot_chap_tot")),
                             tabPanel("Free", align = "center", plotlyOutput("plot_chap_free")),
                             tabPanel("% bound", align = "center", plotlyOutput("plot_chap_bound")),
                             tabPanel("% target and chaperone bound", align = "center", plotlyOutput("plot_Cbnd_Tbnd")),
                             tabPanel("% ternary complex", align = "center", plotlyOutput("plot_PLE_E")),
                             tabPanel("Tracer complex", align = "center", plotlyOutput("plot_EC_E_E_Ctot")), 
                             tabPanel("Tracer", align = "center", plotlyOutput("plot_EC_E"))
                 )
        ),
        tabPanel("Parameters",
                 tableOutput("param")
                 ),
        tabPanel("Info", 
          br(),
          hr(),
          p("This app was developed by ", strong("Richard Stein"), " (EDA-PMX, richard.stein@novartis.com) in collaboration with", strong("Hans-Jörg Roth"), " and team (GDC, hans-joerg.roth@novartis.com) based on Douglass et al. (EF Douglass Jr. et al., 2013: A comprehensive mathematical model for three-body binding equilibria.", em("J. Am. Chem. Soc."), "135(16):6092–6099), and Bomie Han's paper (B Han, 2020: A suite of mathematical solutions to describe ternary complex formation and their application to targeted protein degradation by hetero-bifunctional ligands.", em("J. Biol. Chem."), "295(45):15280–15291)."),
          p("\n"),
          p("\n"),
          div(img(src = "Scheme_Ternary_detailed.png", alt = "Scheme", width = 0.65*1050 , height = 0.65*389), style="text-align: center;"),
          p("\n"),
          p("The formation of ternary complex \\(\\text{[CLT]}\\) through binary complexes can occur along three pathways:"),
          tags$ul(
            tags$li(tags$b("Ligand–target"), "\\(\\text{[LT]}\\) complex formation:"),
            p("\\(K_\\text{CLT} = \\frac{\\text{[L]}\\cdot \\text{[T]}}{\\text{[LT]}} \\cdot \\frac{\\text{[LT]} \\cdot \\text{[C]}}{\\text{[CLT]}} = \\frac{\\text{[C]} \\cdot \\text{[L]} \\cdot \\text{[T]}}{\\text{[CLT]}} = K_\\text{T,1}\\cdot K_\\text{C,2},
              \\)", style="text-align: center;"),
            tags$li(tags$b("Ligand–chaperone"), "\\(\\text{[LC]}\\) complex formation:"),
            p("\\(K_\\text{CLT} = \\frac{\\text{[L]}\\cdot \\text{[C]}}{\\text{[LC]}} \\cdot \\frac{\\text{[LC]} \\cdot \\text{[T]}}{\\text{[CLT]}} = \\frac{\\text{[C]} \\cdot \\text{[L]} \\cdot \\text{[T]}}{\\text{[CLT]}} = K_\\text{C,1}\\cdot K_\\text{T,2},
              \\)", style="text-align: center;"),
            tags$li(tags$b("Target–chaperone"), "\\(\\text{[TC]}\\) complex formation (for intrinsic interactions):"),
            p("\\(K_\\text{CLT} = \\frac{\\text{[T]}\\cdot \\text{[C]}}{\\text{[TC]}} \\cdot \\frac{\\text{[TC]} \\cdot \\text{[L]}}{\\text{[CLT]}} = \\frac{\\text{[C]} \\cdot \\text{[L]} \\cdot \\text{[T]}}{\\text{[CLT]}} = K_\\text{TC}\\cdot K_\\text{L}.
              \\)", style="text-align: center;"),
          ),
          p("Irrespective of the path, all three above described pathways are leading to the formation of the same ternary complex (i.e., \\( K_\\text{C,1} \\cdot K_\\text{T,2} = K_\\text{T,1} \\cdot K_\\text{C,2} = K_\\text{TC} \\cdot K_\\text{L}\\) ). Based on this, we define three", strong("cooperativies"), "describing:
          "),
          tags$ul(
            tags$li("The influence of the counter protein on the two ligand–protein complexes (\\(\\text{[CL]}\\) or \\(\\text{[TL]}\\)):"),
            p("\n"),
            p("\\( \\alpha  = \\frac{K_\\text{T,1}}{K_\\text{T,2}} = \\frac{K_\\text{C,1}}{K_\\text{C,2}}. 
              \\)", style="text-align: center;"),
            p("The cooperativity \\( \\alpha \\) is defined as the ratio of the affinities of free ligand \\(\\text{L}\\) to target \\(\\text{T}\\) in the absence and presence of prebound chaperone to ligand, \\(K_\\text{T,1}\\) and \\(K_\\text{T,2}\\), respectively, or the ratio of the affinities of free ligand \\(\\text{L}\\) to chaperone \\(\\text{C}\\) in the absence and presence of prebound target to ligand, \\(K_\\text{C,1}\\) and \\(K_\\text{C,2}\\), respectively."),
            tags$li("The influence of a preformed protein–protein complex (for intrinsic interactions) on the ligand binding to either of the two proteins:"),
            p("\n"),
            p("\\( \\beta  = \\frac{K_\\text{T,1}}{K_\\text{L}} = \\frac{K_\\text{TC}}{K_\\text{C,2}},\\quad \\gamma  = \\frac{K_\\text{C,1}}{K_\\text{L}} = \\frac{K_\\text{TC}}{K_\\text{T,2}}. 
              \\)", style="text-align: center;"),
            p("Here, the cooperativity \\( \\beta \\) is defined as the ratio of the affinities of free ligand \\(\\text{L}\\) to target \\(\\text{T}\\) in the absence and presence of prebound target to chaperone, \\(K_\\text{T,1}\\) and \\(K_\\text{L}\\), respectively, or the ratio of the affinities of target \\(\\text{T}\\) to chaperone \\(\\text{C}\\) in the absence and presence of prebound ligand to target, \\(K_\\text{TC}\\) and \\(K_\\text{C,2}\\), respectively. Analogously, the cooperativity \\( \\gamma \\) is defined as the ratio of the affinities of free ligand \\(\\text{L}\\) to chaperone \\(\\text{C}\\) in the absence and presence of prebound target to chaperone, \\(K_\\text{C,1}\\) and \\(K_\\text{L}\\), respectively, or the ratio of the affinities of target \\(\\text{T}\\) to chaperone \\(\\text{C}\\) in the absence and presence of prebound ligand to chaperone, \\(K_\\text{TC}\\) and \\(K_\\text{T,2}\\), respectively."),
          ),
          p(strong("Many of the presented results including the source code of R-ShinyApp can be found in our publication:"), "RR Stein, M Fouché, JD Kearns, HJ Roth, 2023: A model-informed method to retrieve intrinsic from apparent cooperativity and project cellular target occupancy for ternary complex-forming compounds.", em("RSC Chem Biol."), "19;4(7):512-523.
          ")
        )
        #,
        #p("© Richard Stein (PKS M&S)")
      )
      ,
      p("\n"),
      p("© Richard Stein (EDA-PMX)")
    )
  )
)
