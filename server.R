# Ternary complex steady state concentration ShinyApp

rm(list=ls())

library(shiny)
library(shinyjs)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
require(deSolve)
require(rootSolve)
library(reshape)
library(plotly)

# compute_PLE <- function(Lfree_Cfree, param) {
#   if (class(Lfree_Cfree) == "matrix" | class(Lfree_Cfree) == "data.frame") {
#     Lfree <- Lfree_Cfree[,1]
#     Cfree <- Lfree_Cfree[,2]
#   } 
#   else {
#     Lfree <- Lfree_Cfree[1]
#     Cfree <- Lfree_Cfree[2]
#   }
#   
#   if (param$Tracer_binding[1] == 'competitive') {
#     KEprime <- ifelse(param$KC > 0, Cfree/param$KC, 0)
#   } 
#   else {
#     # in this case: Cfree = 0
#     KEprime <- param$Kconf
#     #Lfree <- Lfree_Cfree
#   }  
#   
#   f <- param$Ptot + param$Etot + 1/(param$alpha * Lfree) * (Lfree + param$KP1) * (Lfree + (param$KE1*(1+KEprime)))
#   PLE <- ( f - (f^2-4*param$Ptot*param$Etot)^(1/2) )/2   
#   PL <- (param$Ptot - PLE) * Lfree/(Lfree + param$KP1)
#   EL <- (param$Etot - PLE) * Lfree/(Lfree + (param$KE1*(1+KEprime)))
#   P <- (param$Ptot - PLE) * param$KP1/(Lfree + param$KP1)
#   E <- (param$Etot - PLE) * param$KE1/(Lfree + param$KE1*(1+KEprime)) # this is free E! No EC or Ec included
#   
#   return(data.frame(cbind(PLE = PLE, PL = PL, EL = EL, P = P, E = E)))
# }

# Compute Ltot as a function of param and PLE => TI50, TF50
compute_Ltot_from_PLE <- function(PLE, param){
  with(as.list(c(PLE, param)), {
    tmp1 <- param$Ptot + param$Etot - (param$KP1 + param$KE1) + param$alpha*(PLE-param$Ptot)*(PLE-param$Etot)/PLE + param$Ptot*(param$KE1-param$KP1)/param$alpha/(PLE-param$Ptot) + param$Etot*(param$KP1-param$KE1)/param$alpha/(PLE-param$Etot)
    tmp2 <- ( param$alpha/PLE + 1/(param$Ptot-PLE) + 1/(param$Etot-PLE) ) * ( (param$alpha*(PLE-param$Ptot)*(PLE-param$Etot)-PLE*param$KP1)^2 - 2*PLE*param$KE1*(param$alpha*(PLE-param$Ptot)*(PLE-param$Etot)+PLE*param$KP1) + PLE^2*param$KE1^2 )^(0.5)/param$alpha
    
    # Left/- solution to Ltot_l = f(PLE; param) (from Eq. (6) in Douglass et al.)
    Ltot_l <- 0.5*(tmp1 - tmp2) 
    # Right/+ solution to Ltot_r = f(PLE; param) (from Eq. (6) in Douglass et al.)
    Ltot_r <- 0.5*(tmp1 + tmp2) 
    
    return(data.frame(Ltot_l, Ltot_r))
  })
}

# Objective function
#func_obj <- function(Lfree, param){
#  val_tmp <- compute_PLE(Lfree, param)
#  val <- (Lfree + val_tmp$PL + val_tmp$EL + val_tmp$PLE - param$Ltot)^2
#  return(cbind(f = val))
#}

# func_obj <- function(Lfree_Cfree, param){
#   if (param$Tracer_binding[1] == "competitive") {
#     val_tmp <- compute_PLE(Lfree_Cfree, param)
#     Lfree <- Lfree_Cfree[1]
#     Cfree <- Lfree_Cfree[2]
#     val1 <- (Lfree + val_tmp$PL + val_tmp$EL + val_tmp$PLE - param$Ltot)^2
#     val2 <- (Cfree + ifelse(param$KC > 0, Cfree*val_tmp$E/param$KC, 0) - param$Ctot)^2
#     val = c(val1, val2)
#   } else {
#     Lfree <- Lfree_Cfree
#     val_tmp <- compute_PLE(Lfree, param)
#     val <- (Lfree + val_tmp$PL + val_tmp$EL + val_tmp$PLE - param$Ltot)^2
#   }
#   return(cbind(f = val))
# }

Find_root <- function(param){
  Efree_Pfree_Lfree_min <- matrix(0,length(param$Ltot),3)
  
  for (id in 1:length(param$Ltot)) {
    Efree0 <- param$Etot[1]
    Pfree0 <- param$Ptot[1]
    Lfree0 <- param$Ltot[id]
    Efree_Pfree_Lfree0 <- cbind(Efree0, Pfree0, Lfree0)
    
    root <- multiroot(f = func_obj_intrinsic, start = Efree_Pfree_Lfree0, parms = param[id,], verbose = FALSE, useFortran = TRUE, positive = TRUE)$root
    Efree_Pfree_Lfree_min[id,] <- c(root[1], root[2], root[3])
  }
  
  res <- param %>% mutate(Efree = Efree_Pfree_Lfree_min[,1], Pfree = Efree_Pfree_Lfree_min[,2], Lfree = Efree_Pfree_Lfree_min[,3]) %>% 
    mutate(PE = Efree*Pfree * ifelse(param$KPE > 0, 1/param$KPE, 0), # Efree*Pfree/param$KPE, 
           EL = Efree*Lfree * ifelse(param$KE1 > 0, 1/param$KE1, 0), 
           EC = 0,
           PL = Pfree*Lfree * ifelse(param$KP1 > 0, 1/param$KP1, 0), 
           PLE = Pfree*Lfree*Efree * ifelse(param$KP1 > 0, 1/param$KP1, 0) * ifelse(param$KE1 > 0, 1/param$KE1, 0) * param$alpha) 
  
  #res <- res[rowSums(is.na(res)) == 0,]
  return(f = res)
}

func_obj_intrinsic <- function(Efree_Pfree_Lfree, param){
  Efree <- Efree_Pfree_Lfree[1]
  Pfree <- Efree_Pfree_Lfree[2]
  Lfree <- Efree_Pfree_Lfree[3]
  
  PE <- ifelse(param$KPE > 0, Efree*Pfree/param$KPE, 0) # /param$KPE
  EL <- ifelse(param$KE1 > 0, Efree*Lfree/param$KE1, 0) # /param$KE1
  PL <- ifelse(param$KP1 > 0, Pfree*Lfree/param$KP1, 0) # /param$KP1
  PLE <- ifelse(param$KP1 > 0, Pfree*Lfree*Efree/param$KP1, 0) * ifelse(param$KE1 > 0, 1/param$KE1, 0) * param$alpha # /param$KP1/param$KE1
  
  val1 <- (Efree+EL+PE+PLE-param$Etot)^2 # difference of Etot
  val2 <- (Pfree+PL+PE+PLE-param$Ptot)^2 # difference of Ptot
  val3 <- (Lfree+EL+PL+PLE-param$Ltot)^2 # difference of Ltot
  val = c(val1, val2, val3)
  
  return(f = val)
}

func_obj_intrinsic_tracer <- function(Efree_Pfree_Lfree_Cfree, param){
  Efree <- Efree_Pfree_Lfree_Cfree[1]
  Pfree <- Efree_Pfree_Lfree_Cfree[2]
  Lfree <- Efree_Pfree_Lfree_Cfree[3]
  Cfree <- Efree_Pfree_Lfree_Cfree[4]
    
  PE <- ifelse(param$KPE > 0, Efree*Pfree/param$KPE, 0) # /param$KPE
  EL <- ifelse(param$KE1 > 0, Efree*Lfree/param$KE1, 0) # /param$KE1
  PL <- ifelse(param$KP1 > 0, Pfree*Lfree/param$KP1, 0) # /param$KP1
  EC <- ifelse(param$KC > 0, Cfree*Efree/param$KC, 0)
  PLE <- ifelse(param$KP1 > 0, Pfree*Lfree*Efree/param$KP1, 0) * ifelse(param$KE1 > 0, 1/param$KE1, 0) * param$alpha # /param$KP1/param$KE1
  
  val1 <- (Efree+EL+PE+EC+PLE-param$Etot)^2 # difference of Etot
  val2 <- (Pfree+PL+PE+PLE-param$Ptot)^2 # difference of Ptot
  val3 <- (Lfree+EL+PL+PLE-param$Ltot)^2 # difference of Ltot
  val4 <- (Cfree+EC-param$Ctot)^2
  val = c(val1, val2, val3, val4)
  
  return(f = val)
}

funct_obj_CboundEC50 <- function(Ltot, param){
  param$Ltot <- Ltot
  res <- Find_root(param[1,])
  # Cbound at EC50 is supposed to be 50%!
  res %>% mutate(Ebound = ifelse(param$Etot>0, (EL+PE+ifelse("EC" %in% names(param), EC, 0)+PLE)/param$Etot*100, 0)) %>% 
    select(Ebound) %>% unlist %>% as.numeric -> Cbound_EC50_est
  diff <- (Cbound_EC50_est - 50)^2
  
  return(f = diff)
}

funct_obj_TboundEC50 <- function(Ltot, param){
  param$Ltot <- Ltot
  res <- Find_root(param[1,])
  # Tbound at EC50 is supposed to be 50%!
  res %>% mutate(Pbound = ifelse(param$Ptot>0, (PL+PE+PLE)/param$Ptot*100, 0)) %>% 
    select(Pbound) %>% unlist %>% as.numeric -> Tbound_EC50_est
  diff <- (Tbound_EC50_est - 50)^2
  
  return(f = diff)
}

Find_EC50_Cbound_Tbound <- function(param){
  Ltot_Cbound_EC50 <- nlm(funct_obj_CboundEC50, param[1,]$KE1+param[1,]$Etot/2, param[1,])$estimate
  Ltot_Tbound_EC50 <- nlm(funct_obj_TboundEC50, param[1,]$KP1+param[1,]$Ptot/2, param[1,])$estimate
  
  return(as.numeric(c(Ltot_Cbound_EC50 = Ltot_Cbound_EC50, Ltot_Tbound_EC50 = Ltot_Tbound_EC50)))
}

server <- function(input, output, session) {

  observeEvent(input$Native_interaction, {
    if(input$Native_interaction == 'none') {
      #useShinyjs()  # Set up shinyjs
      shinyjs::enable("alpha")
      updateNumericInput(session, "KPE", value = 0)
      shinyjs::disable("KPE")
      updateNumericInput(session, "KL", value = 0)
      if (!is.na(input$Cooperativity) && input$Cooperativity_choice == "alpha") {
        updateNumericInput(session, "alpha", value = input$Cooperativity)
      } else {
        shinyjs::enable("alpha")
        shinyjs::disable("Cooperativity")
      }
      updateNumericInput(session, "Cooperativity", value = NA)
      shinyjs::disable("KPE")
      shinyjs::disable("Cooperativity")
    } 
    else if(input$Native_interaction == 'native') {
        #shinyjs::useShinyjs()  # Set up shinyjs
        shinyjs::enable("KPE")
        shinyjs::enable("Cooperativity")
        shinyjs::disable("alpha")
        updateNumericInput(session, "KPE", value = 1000)
        updateNumericInput(session, "Cooperativity", value = input$alpha)
        updateNumericInput(session, "alpha", value = NA)
        shinyjs::disable("alpha")
    } 
    #else {
    #  shinyjs::useShinyjs()  # Set up shinyjs
    #  shinyjs::enable("alpha")
    #  updateNumericInput(session, "KL", value = 0)
    #  updateNumericInput(session, "KPE", value = 0)
    #  updateNumericInput(session, "Cooperativity", value = NA)
    #  shinyjs::disable("Cooperativity", "KPE")
    #}    
  })
  
  observeEvent(input$Tracer_binding, {
    if(input$Tracer_binding == 'none'){
      updateNumericInput(session, "KC", value = 0)
      updateNumericInput(session, "Ctot", value = 0)
      updateNumericInput(session, "Kconf", value = 0)
    }else if(input$Tracer_binding == 'competitive'){
      updateNumericInput(session, "KC", value = 65)
      updateNumericInput(session, "Ctot", value = 25)
      updateNumericInput(session, "Kconf", value = 0)
    }else if(input$Tracer_binding == 'conf_eq'){
      updateNumericInput(session, "KC", value = 0)
      updateNumericInput(session, "Ctot", value = 0)
      updateNumericInput(session, "Kconf", value = 100)
    }
  })
    
  res <- reactive({
    
    validate(
      need(as.numeric(unlist(strsplit(input$Ligtot, "\n"))) + as.numeric(unlist(strsplit(input$Tartot, "\n"))) + as.numeric(unlist(strsplit(input$E3tot, "\n")))> 0, "Please provide ligand, target and E3/FKBP concentration(s)!"),
      #need(input$alpha > 0, "Please provide α>0!"),
      need(input$KP1 >= 0, "Please provide KT1!"),
      need(input$KE1 >= 0, "Please provide KC1!"),
      need(input$KC >= 0, "Please provide KC!"),
      need(input$Kconf >= 0, "Please provide Kconf!"),
      #need(input$KPE >= 0, "Please provide KTC!"),
      #need(input$Cooperativity > 0, "Please provide Cooperativity!")
      if(input$Native_interaction == 'none') {
        need(input$alpha > 0, "Please provide α>0!")
      } else if(input$Native_interaction == 'native' && input$KT_choice == 'KT2' && input$Cooperativity_choice == 'gamma') 
        need(input$alpha > 0 & input$Cooperativity>0 & input$KPE>0, "Please provide α>0, γ>0 and/or KTC > 0!")
      else
      {
        need(input$Cooperativity>0 & input$KPE>0, "Please provide Cooperativity and/or KTC > 0!")
      }
    )
    #req(input$Ligtot)
    #req(input$Tartot) 
    #req(input$alpha) 
    
    if(input$Cooperativity_choice == 'alpha' && input$Native_interaction == 'native') {
      alpha_comp = input$Cooperativity
    } else if(input$Cooperativity_choice == 'beta' && input$Native_interaction == 'native') {
      #KL = input$KP1/input$Cooperativity
      alpha_comp = input$KE1 * input$Cooperativity/input$KPE # input$KP1*input$KE1/(input$KP1/input$Cooperativity)/input$KPE
    } else if(input$Cooperativity_choice == 'gamma' && input$Native_interaction == 'native') {
      #KL = input$KC1/input$Cooperativity
      # does not work for input$KT_choice == "KT2"!
      alpha_comp = input$KP1*input$Cooperativity/input$KPE # input$KP1*input$KE1/(input$KE1/input$Cooperativity)/input$KPE
    } else {
      alpha_comp = input$alpha
    }
    
    if(input$KT_choice == 'KT2') {
      KP1_comp = input$KP1*alpha_comp
    }
    else {
      KP1_comp = input$KP1
    }
    
    param = data.frame(Ltot = as.numeric(unlist(strsplit(input$Ligtot, "\n"))), # Total ligand 
                       Ptot = as.numeric(unlist(strsplit(input$Tartot, "\n"))), # Total target
                       Etot = as.numeric(unlist(strsplit(input$E3tot, "\n"))), # Total chaperone
                       KP1 = KP1_comp, # input$KP1, # KP is in nM
                       KE1 = input$KE1,
                       KPE = input$KPE,
                       #KL = input$KL,
                       Kconf = input$Kconf, 
                       Ctot = input$Ctot,
                       KC = input$KC,
                       alpha = alpha_comp,
                       #alpha = input$alpha,
                       Tracer_binding = input$Tracer_binding,
                       Native_interaction = input$Native_interaction
                       )

    # Add EC50s to Ltot 
    Ltot_EC50_vec <- Find_EC50_Cbound_Tbound(param)
    
    EC50_vec <- Find_EC50_Cbound_Tbound(param)
    param <- rbind(param, 
                      data.frame(Ltot = as.numeric(Ltot_EC50_vec[1]), param[1,2:length(param)]),
                      data.frame(Ltot = as.numeric(Ltot_EC50_vec[2]), param[1,2:length(param)]))
    
    # Sort param
    param <- arrange(param, Ltot, Ptot, Etot)
    
    # Root finding: Use non-linear root finder for all cases (with or without native interaction)
    # no tracer
    if(param$Tracer_binding[1] == 'none') {
      res <- Find_root(param)
      # Efree_Pfree_Lfree_min <- matrix(0,length(param$Ltot),3)
      # 
      # for (id in 1:length(param$Ltot)) {
      #   Efree0 <- param$Etot[1]
      #   Pfree0 <- param$Ptot[1]
      #   Lfree0 <- param$Ltot[id]
      #   Efree_Pfree_Lfree0 <- cbind(Efree0, Pfree0, Lfree0)
      #   
      #   root <- multiroot(f = func_obj_intrinsic, start = Efree_Pfree_Lfree0, parms = param[id,], verbose = FALSE, useFortran = TRUE, positive = TRUE)$root
      #   Efree_Pfree_Lfree_min[id,] <- c(root[1], root[2], root[3])
      # }
      # 
      # res <- param %>% mutate(Efree = Efree_Pfree_Lfree_min[,1], Pfree = Efree_Pfree_Lfree_min[,2], Lfree = Efree_Pfree_Lfree_min[,3], Cfree = 0) %>% 
      #   mutate(PE = ifelse(param$KPE > 0, Efree*Pfree/param$KPE, 0), # Efree*Pfree/param$KPE, 
      #          EL = ifelse(param$KE1 > 0, Efree*Lfree/param$KE1, 0),
      #          EC = 0,
      #          PL = ifelse(param$KP1 > 0, Pfree*Lfree/param$KP1, 0), 
      #          PLE = ifelse(param$KP1 > 0, Pfree*Lfree*Efree/param$KP1, 0) * ifelse(param$KE1 > 0, 1/param$KE1, 0) * param$alpha
      #          ) 

    } 
    # with tracer or conformational equilibrium (the latter to be implemented!)
    else {
      Efree_Pfree_Lfree_Cfree_min <- matrix(0,length(param$Ltot),4)
      
      for (id in 1:length(param$Ltot)) {
        Efree0 <- param$Etot[1]
        Pfree0 <- param$Ptot[1]
        Lfree0 <- param$Ltot[id]
        Cfree0 <- param$Ctot[1]
        Efree_Pfree_Lfree_Cfree0 <- cbind(Efree0, Pfree0, Lfree0, Cfree0)
        
        root <- multiroot(f = func_obj_intrinsic_tracer, start = Efree_Pfree_Lfree_Cfree0, parms = param[id,], verbose = FALSE, useFortran = TRUE, positive = TRUE)$root
        Efree_Pfree_Lfree_Cfree_min[id,] <- c(root[1], root[2], root[3], root[4])
      }
      
      res <- param %>% mutate(Efree = Efree_Pfree_Lfree_Cfree_min[,1], Pfree = Efree_Pfree_Lfree_Cfree_min[,2], Lfree = Efree_Pfree_Lfree_Cfree_min[,3], Cfree = Efree_Pfree_Lfree_Cfree_min[,4]) %>% 
        # PE: TC target--chaperone, EL: CL chaperone--ligand, EC: chaperone--tracer, PL: TL target--ligand, PLE: CLT ternary complex 
        mutate(PE = ifelse(param$KPE > 0, Efree*Pfree/param$KPE, 0), # Efree*Pfree/param$KPE, 
               EL = ifelse(param$KE1 > 0, Efree*Lfree/param$KE1, 0), 
               EC = ifelse(param$KC > 0, Efree*Cfree/param$KC, 0),
               PL = ifelse(param$KP1 > 0, Pfree*Lfree/param$KP1, 0), 
               PLE = ifelse(param$KP1 > 0, Pfree*Lfree*Efree/param$KP1, 0) * ifelse(param$KE1 > 0, 1/param$KE1, 0) * param$alpha
        )
    }
    
# remove rows with NA
    #res[is.na(res)] <- 0   
    res <- res[rowSums(is.na(res)) == 0,]
    
    # Maximum ternary complex concentration
    PLEmax <- 0.5*(param$Ptot + param$Etot + ((param$KP1)^(0.5) + (param$KE1)^(0.5))^2/param$alpha - ( (param$Ptot + param$Etot + ((param$KP1)^(0.5) + (param$KE1)^(0.5))^2/param$alpha)^2 - 4*param$Ptot*param$Etot )^(0.5))[1]
    
    # TF50 and TI50
    T50 <- compute_Ltot_from_PLE(PLEmax/2, param[1,])
    Ltotmax <- (1/((param$KP1)^(0.5) + (param$KE1)^(0.5)) * ((param$KE1)^(0.5)*(param$Ptot+param$KP1)+(param$KP1)^(0.5)*(param$Etot+param$KE1)))[1]
    Lfreemax <- ((param$KP1*param$KE1)^(0.5))[1]
    
    # EC_Lsmall
    EC_Lsmall <- (param$Ctot+param$Etot+param$KC)/2 - ( (param$Ctot+param$Etot+param$KC)^2/4 - param$Ctot*param$Etot)^(0.5)    
    
    # Ligand output 
    if(input$Ligand_display == 'Free_ligand') {
      Display_Ltot = 0
    } else if (input$Ligand_display == 'Total_ligand') {
      Display_Ltot = 1
    } else if (input$Ligand_display == 'Total_target') {
      Display_Ltot = 2
    } else if (input$Ligand_display == 'Total_chap') {
      Display_Ltot = 3
    }
        
    # Combine all data for plotting
    res <- res %>% mutate(PLEmax, 
                          TF50 = as.numeric(T50[1]), 
                          TI50 = as.numeric(T50[2]), 
                          Ltotmax = as.numeric(Ltotmax),
                          Lfreemax = ifelse(!is.na(as.numeric(Lfreemax)), as.numeric(Lfreemax), 0),
                          EC_Lsmall = as.numeric(EC_Lsmall)[1],
                          # Ligand
                          L_bnd = (PL+EL+PLE)/param$Ltot*100,
                          # Target
                          P_bnd = ifelse(param$Ptot>0, round((PL+PE+PLE)/param$Ptot*100,2), 0), # TO
                          PLE_Ptot = ifelse(param$Ptot>0, PLE/param$Ptot*100, 0),
                          Ltot_EC50_Tbound = as.numeric(EC50_vec[2]),
                          # Chaperone
                          E_bnd = ifelse(param$Etot>0, round((EL+PE+EC+PLE)/param$Etot*100,2), 0),
                          PLE_Etot = ifelse(param$Etot>0, PLE/param$Etot*100, 0),
                          Ltot_EC50_Cbound = as.numeric(EC50_vec[1]),
                          # Display
                          Display_Ltot = as.numeric(Display_Ltot[1]),
                          Display_PLEmax = as.numeric(ifelse(param$Tracer_binding[1] == "none" & param$Native_interaction[1] == "none", 1, 0)), 
                          )
    
    return(res)
  })

  # Plot output
  #group.colors <- c("PLE" = "#F8766D", "PL" = "#00BA38", "EL" ="#619CFF", "Ltot" = "#C77CFF", "Lfree" = "#C77CFF", "L_bnd" = "#C77CFF", "Ptot" = "#00BA38", "P" = "#00BA38", "P_bnd" = "#00BA38", "PLE_Ptot" = "#F8766D", "Etot" = "#619CFF", "E" = "#619CFF", "E_bnd" = "#619CFF", "PLE_Etot" = "#F8766D", "Cfree" = "#00BFC4", "EC" = "#FB61D7", "Ctot" = "#A58AFF")
  group.colors <- c("TLC" = "#F8766D", "TL" = "#00BA38", "TC" = "#C4961A", "CL" = "#619CFF", "Ltot" = "#C77CFF", "Lfree" = "#C77CFF", "L_bound" = "#C77CFF", "Ttot" = "#00BA38", "Tfree" = "#00BA38", "T_bound" = "#00BA38", "TLC_Ttot" = "#F8766D", "Ctot" = "#619CFF", "Cfree" = "#619CFF", "C_bound" = "#619CFF", "TLC_Ctot" = "#F8766D", "Trfree" = "#00BFC4", "CTr" = "#FB61D7", "Trtot" = "#A58AFF")
  
  # Complex combined
  
  # Only ggplot (no plotly)
  # output$plot_cmplx <- renderPlot({
  #   res = res()
  # 
  #   if (res$Display_Ltot[1] == 1) {
  #     pk <- res %>% select(Ltot, PLE, PL, EL) %>% melt(id.vars = 'Ltot')  %>% dplyr::rename(xval = Ltot)
  #   } else if (res$Display_Ltot[1] == 0) {
  #     pk <- res %>% select(Lfree, PLE, PL, EL) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
  #   } else if (res$Display_Ltot[1] == 2) {
  #     pk <- res %>% select(Ptot, PLE, PL, EL) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot)
  #   } else if (res$Display_Ltot[1] == 3) {
  #     pk <- res %>% select(Etot, PLE, PL, EL) %>%  melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot)
  #   }
  # 
  #   pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable))  + geom_line(alpha=1, size=1.25) +
  #     #annotation_logticks(sides = "b") + #scale_x_continuous(trans="log10", #limits=c(1e-1,NA),
  #                                       #                    breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                                       #                    labels = scales::trans_format("log10", scales::math_format(10^.x))
  #     #) +
  #     scale_color_manual(name = NULL, values = group.colors, labels = c("Ternary complex PLE", "Target–Ligand PL", "Chaperone–Ligand EL")) 
  #     #+ theme(text = element_text(size=14))
  # 
  #   if (res$Display_Ltot[1] == 1){
  #     pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Concentration of complex (in nM)")
  #   } else if (res$Display_Ltot[1] == 0) {
  #     pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Concentration of complex (in nM)")
  #   } else if (res$Display_Ltot[1] == 2) {
  #     pk <- pk + labs(x = "Total target concentration (in nM)", y = "Concentration of complex (in nM)")
  #   } else if (res$Display_Ltot[1] == 3) {
  #     pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Concentration of complex (in nM)")
  #   }
  # 
  #   print(pk)
  # })
  
  # With plotly
  output$plot_cmplx <- renderPlotly({
    res = res()

    if (res$Display_Ltot[1] == 1) {
      pk <- res %>% select(Ltot, PLE, PL, EL, PE) %>% dplyr::rename(TLC = PLE, TL = PL, CL = EL, TC = PE) %>% melt(id.vars = 'Ltot') %>% dplyr::rename(xval = Ltot)
    } else if (res$Display_Ltot[1] == 0) {
      pk <- res %>% select(Lfree, PLE, PL, EL, PE) %>% dplyr::rename(TLC = PLE, TL = PL, CL = EL, TC = PE) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
    } else if (res$Display_Ltot[1] == 2) {
      pk <- res %>% select(Ptot, PLE, PL, EL, PE) %>% dplyr::rename(TLC = PLE, TL = PL, CL = EL, TC = PE) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot)
    } else if (res$Display_Ltot[1] == 3) {
      pk <- res %>% select(Etot, PLE, PL, EL, PE) %>% dplyr::rename(TLC = PLE, TL = PL, CL = EL, TC = PE) %>% melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot)
    }
    
    pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) + 
      scale_x_continuous(trans="log10", breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      scale_color_manual(name = NULL, labels = c(1,2,3), values = group.colors)
    
    if (res$Display_Ltot[1] == 1){
      pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Concentration of complex (in nM)")
    } else if (res$Display_Ltot[1] == 0) {
      pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Concentration of complex (in nM)")
    } else if (res$Display_Ltot[1] == 2) {
      pk <- pk + labs(x = "Total target concentration (in nM)", y = "Concentration of complex (in nM)")
    } else if (res$Display_Ltot[1] == 3) {
      pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Concentration of complex (in nM)")
    }
    
    pk <- ggplotly(pk) # %>% layout(xaxis = list(exponentformat = "power"))
    print(pk)
  })
  
  
  # Ternary complex
  
  # Only ggplot (no plotly)
  # output$plot_PLE <- renderPlot({
  #   res = res()
  #   PLEmax <- res$PLEmax[1]
  #   TF50 <- res$TF50[1]
  #   TI50 <- res$TI50[1]
  #   Ltotmax <- res$Ltotmax[1]
  #   Lfreemax <- res$Lfreemax[1]
  # 
  #   if (res$Display_Ltot[1] == 1) { 
  #     pk <- res %>% select(Ltot, PLE) %>% melt(id.vars = 'Ltot')  %>% dplyr::rename(xval = Ltot)
  #   } else if (res$Display_Ltot[1] == 0) { 
  #     pk <- res %>% select(Lfree, PLE) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
  #   } else if (res$Display_Ltot[1] == 2) { 
  #     pk <- res %>% select(Ptot, PLE) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot)
  #   } else if (res$Display_Ltot[1] == 3) { 
  #     pk <- res %>% select(Etot, PLE) %>%  melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot)
  #   }    
  #   
  #   pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) +
  #       annotation_logticks(sides = "b") + scale_x_continuous(trans="log10", #limits=c(1e-1,NA),
  #                                                           breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                                                           labels = scales::trans_format("log10", scales::math_format(10^.x))
  #     ) + theme(text = element_text(size=14)) +
  #     scale_color_manual(name = NULL, values = group.colors)
  #   
  #   if (res$Display_Ltot[1] == 1){ 
  #     pk <- pk + geom_hline(yintercept = as.numeric(PLEmax), linetype = "dashed", size = .75) + annotate("text", x = res$Ltot %>% last(), y = as.numeric(PLEmax), hjust = 0.75, vjust = 1.25, label = "PLE[max]", size = 6, parse = TRUE) + 
  #       geom_vline(xintercept = as.numeric(TF50), linetype = "dashed", size = .75) + annotate("text", x = as.numeric(TF50), y = res$PLE %>% min(), hjust = 0, vjust = -0.25, angle = 90, label = "TF[50]", size = 6, parse = TRUE) + 
  #       geom_vline(xintercept = as.numeric(TI50), linetype = "dashed", size = .75) + annotate("text", x = as.numeric(TI50), y = res$PLE %>% min(), hjust = -0.05, vjust = 1.25, angle = 90, label = "TI[50]", size = 6, parse = TRUE) + 
  #       geom_vline(xintercept = as.numeric(Ltotmax), linetype = "dashed", size = .75) + annotate("text", x = as.numeric(Ltotmax), y = res$PLE %>% min(), hjust = -0.1, vjust = -0.25, angle = 90, label = "L['tot,max']", size = 6, parse = TRUE) + 
  #       labs(x = "Total ligand concentration (in nM)", y = "Concentration of ternary complex (in nM)")
  #   } else if (res$Display_Ltot[1] == 0) { 
  #     pk <- pk + geom_hline(yintercept = as.numeric(PLEmax), linetype = "dashed", size = .75) + annotate("text", x = res$Ltot %>% last(), y = as.numeric(PLEmax), hjust = 0.75, vjust = 1.25, label = "PLE[max]", size = 6, parse = TRUE) + 
  #       geom_vline(xintercept = as.numeric(Lfreemax), linetype = "dashed", size = .75) + annotate("text", x = as.numeric(Lfreemax), y = res$PLE %>% min(), hjust = -0.1, vjust = -0.25, angle = 90, label = "L['max']", size = 6, parse = TRUE) + 
  #       labs(x = "Free ligand concentration (estimated) (in nM)", y = "Concentration of ternary complex (in nM)") 
  #   } else if (res$Display_Ltot[1] == 2) { 
  #     pk <- pk + labs(x = "Total target concentration (in nM)", y = "Concentration of ternary complex (in nM)") 
  #   } else if (res$Display_Ltot[1] == 3) { 
  #     pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Concentration of ternary complex (in nM)") 
  #   }
  #   
  #   print(pk)
  # })
  
  # With plotly
  output$plot_PLE <- renderPlotly({
    res = res()
    PLEmax <- res$PLEmax[1]
    TF50 <- res$TF50[1]
    TI50 <- res$TI50[1]
    Ltotmax <- res$Ltotmax[1]
    Lfreemax <- res$Lfreemax[1]
    
    if (res$Display_Ltot[1] == 1) { 
      pk <- res %>% select(Ltot, PLE) %>% dplyr::rename(TLC = PLE) %>% melt(id.vars = 'Ltot')  %>% dplyr::rename(xval = Ltot)
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- res %>% select(Lfree, PLE) %>% dplyr::rename(TLC = PLE) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- res %>% select(Ptot, PLE) %>% dplyr::rename(TLC = PLE) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot)
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- res %>% select(Etot, PLE) %>% dplyr::rename(TLC = PLE) %>% melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot)
    }    
    
    pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) + 
      scale_x_continuous(trans="log10", breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      scale_color_manual(name = NULL, values = group.colors) 
    
    if (res$Display_Ltot[1] == 1 & res$Display_PLEmax[1] == 1){ 
      pk <- pk + geom_hline(yintercept = as.numeric(PLEmax), linetype = "dashed", size = .75) + annotate("text", x = res$Ltot %>% last()*0.5, y = as.numeric(PLEmax)*0.95, hjust = 0.75, vjust = 1.25, label = "TLCmax", size = 4, parse = TRUE) + 
        geom_vline(xintercept = as.numeric(TF50), linetype = "dashed", size = .75) + annotate("text", x = as.numeric(TF50), y = res$PLE %>% min(), hjust = 0, vjust = -0.25, angle = 90, label = "TF50", size = 4, parse = TRUE) + 
        geom_vline(xintercept = as.numeric(TI50), linetype = "dashed", size = .75) + annotate("text", x = as.numeric(TI50), y = res$PLE %>% min(), hjust = -0.05, vjust = 1.25, angle = 90, label = "TI50", size = 4, parse = TRUE) + 
        geom_vline(xintercept = as.numeric(Ltotmax), linetype = "dashed", size = .75) + annotate("text", x = as.numeric(Ltotmax), y = res$PLE %>% min(), hjust = -0.1, vjust = -0.25, angle = 90, label = "Ltot,max", size = 4, parse = TRUE) + 
        labs(x = "Total ligand concentration (in nM)", y = "Concentration of ternary complex (in nM)")
    } else if (res$Display_Ltot[1] == 1 & res$Display_PLEmax[1] == 0){ 
        pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Concentration of ternary complex (in nM)")
    } else if (res$Display_Ltot[1] == 0 & res$Display_PLEmax[1] == 1) { 
      pk <- pk + geom_hline(yintercept = as.numeric(PLEmax), linetype = "dashed", size = .75) + annotate("text", x = res$Ltot %>% last()*0.5, y = as.numeric(PLEmax)*0.95, hjust = 0.75, vjust = 1.25, label = "PLEmax", size = 4, parse = TRUE) + 
        geom_vline(xintercept = as.numeric(Lfreemax), linetype = "dashed", size = .75) + annotate("text", x = as.numeric(Lfreemax), y = res$PLE %>% min(), hjust = -0.1, vjust = -0.25, angle = 90, label = "Lmax", size = 4, parse = TRUE) + 
        labs(x = "Free ligand concentration (estimated) (in nM)", y = "Concentration of ternary complex (in nM)") 
    } else if (res$Display_Ltot[1] == 0 & res$Display_PLEmax[1] == 0) { 
      pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Concentration of ternary complex (in nM)") 
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- pk + labs(x = "Total target concentration (in nM)", y = "Concentration of ternary complex (in nM)") 
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Concentration of ternary complex (in nM)") 
    }
    
    pk <- ggplotly(pk) # %>% layout(xaxis = list(exponentformat = "power"))
    print(pk)
  })
  
  
  # EL
  
  # Only ggplot (no plotly)
  # output$plot_EL <- renderPlot({
  #   res = res()
  #   
  #   if (res$Display_Ltot[1] == 1) { 
  #     pk <- res %>% select(Ltot, EL) %>% melt(id.vars = 'Ltot')  %>% dplyr::rename(xval = Ltot)
  #   } else if (res$Display_Ltot[1] == 0) { 
  #     pk <- res %>% select(Lfree, EL) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
  #   } else if (res$Display_Ltot[1] == 2) { 
  #     pk <- res %>% select(Ptot, EL) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot)
  #   } else if (res$Display_Ltot[1] == 3) { 
  #     pk <- res %>% select(Etot, EL) %>%  melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot)
  #   }    
  #   
  #   pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) +
  #       annotation_logticks(sides = "b") + scale_x_continuous(trans="log10", #limits=c(1e-1,NA),
  #                                                           breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                                                           labels = scales::trans_format("log10", scales::math_format(10^.x))
  #       ) +
  #       theme(text = element_text(size=14)) +
  #       scale_color_manual(name = NULL, values = group.colors)
  #   
  #   if (res$Display_Ltot[1] == 1){ 
  #     pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Concentration of chaperone–ligand complex (in nM)")
  #   } else if (res$Display_Ltot[1] == 0) { 
  #     pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Concentration of chaperone–ligand complex (in nM)") 
  #   } else if (res$Display_Ltot[1] == 2) { 
  #     pk <- pk + labs(x = "Total target concentration (in nM)", y = "Concentration of chaperone–ligand complex (in nM)") 
  #   } else if (res$Display_Ltot[1] == 3) { 
  #     pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Concentration of chaperone–ligand complex (in nM)") 
  #   }
  #   
  #   print(pk)
  # })
  
  # With plotly
  output$plot_EL <- renderPlotly({
    res = res()
    
    if (res$Display_Ltot[1] == 1) { 
      pk <- res %>% select(Ltot, EL) %>% dplyr::rename(CL = EL) %>% melt(id.vars = 'Ltot')  %>% dplyr::rename(xval = Ltot)
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- res %>% select(Lfree, EL) %>% dplyr::rename(CL = EL) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- res %>% select(Ptot, EL) %>% dplyr::rename(CL = EL) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot)
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- res %>% select(Etot, EL) %>% dplyr::rename(CL = EL) %>%  melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot)
    }    
    
    pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) + 
      scale_x_continuous(trans="log10", breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      scale_color_manual(name = NULL, values = group.colors) 
    
    if (res$Display_Ltot[1] == 1){ 
      pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Concentration of chaperone–ligand complex (in nM)")
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Concentration of chaperone–ligand complex (in nM)") 
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- pk + labs(x = "Total target concentration (in nM)", y = "Concentration of chaperone–ligand complex (in nM)") 
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Concentration of chaperone–ligand complex (in nM)") 
    }
    
    pk <- ggplotly(pk) # %>% layout(xaxis = list(exponentformat = "power"))
    print(pk)
  })
  
  
  # PL
  
  # Only ggplot (no plotly)
  # output$plot_PL <- renderPlot({
  #   res = res()
  #   
  #   if (res$Display_Ltot[1] == 1) { 
  #     pk <- res %>% select(Ltot, PL) %>% melt(id.vars = 'Ltot')  %>% dplyr::rename(xval = Ltot)
  #   } else if (res$Display_Ltot[1] == 0) { 
  #     pk <- res %>% select(Lfree, PL) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
  #   } else if (res$Display_Ltot[1] == 2) { 
  #     pk <- res %>% select(Ptot, PL) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot)
  #   } else if (res$Display_Ltot[1] == 3) { 
  #     pk <- res %>% select(Etot, PL) %>%  melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot)
  #   }   
  #   
  #   pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) +
  #       annotation_logticks(sides = "b") + scale_x_continuous(trans="log10", #limits=c(1e-1,NA),
  #                                                           breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                                                           labels = scales::trans_format("log10", scales::math_format(10^.x))
  #       ) + theme(text = element_text(size=14)) + scale_color_manual(name = NULL, values = group.colors)
  # 
  #   if (res$Display_Ltot[1] == 1){ 
  #     pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Concentration of target–ligand complex (in nM)")
  #   } else if (res$Display_Ltot[1] == 0) { 
  #     pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Concentration of target–ligand complex (in nM)") 
  #   } else if (res$Display_Ltot[1] == 2) { 
  #     pk <- pk + labs(x = "Total target concentration (in nM)", y = "Concentration of target–ligand complex (in nM)") 
  #   } else if (res$Display_Ltot[1] == 3) { 
  #     pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Concentration of target–ligand complex (in nM)") 
  #   }
  #   
  #   print(pk)
  # })
  
  # With plotly
  output$plot_PL <- renderPlotly({
    res = res()
    
    if (res$Display_Ltot[1] == 1) { 
      pk <- res %>% select(Ltot, PL) %>% dplyr::rename(TL = PL) %>% melt(id.vars = 'Ltot')  %>% dplyr::rename(xval = Ltot)
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- res %>% select(Lfree, PL) %>% dplyr::rename(TL = PL) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- res %>% select(Ptot, PL) %>% dplyr::rename(TL = PL) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot)
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- res %>% select(Etot, PL) %>% dplyr::rename(TL = PL) %>%  melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot)
    }   
    
    pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) + 
      scale_x_continuous(trans="log10", breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      scale_color_manual(name = NULL, values = group.colors) 
    
    if (res$Display_Ltot[1] == 1){ 
      pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Concentration of target–ligand complex (in nM)")
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Concentration of target–ligand complex (in nM)") 
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- pk + labs(x = "Total target concentration (in nM)", y = "Concentration of target–ligand complex (in nM)") 
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Concentration of target–ligand complex (in nM)") 
    }
    
    pk <- ggplotly(pk) # %>% layout(xaxis = list(exponentformat = "power"))
    print(pk)
  })
  
  # PE
  
  # Only with plotly!
  output$plot_PE <- renderPlotly({
    res = res()
    
    if (res$Display_Ltot[1] == 1) { 
      pk <- res %>% select(Ltot, PE) %>% dplyr::rename(TC = PE) %>% melt(id.vars = 'Ltot')  %>% dplyr::rename(xval = Ltot)
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- res %>% select(Lfree, PE) %>% dplyr::rename(TC = PE) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- res %>% select(Ptot, PE) %>% dplyr::rename(TC = PE) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot)
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- res %>% select(Etot, PE) %>% dplyr::rename(TC = PE) %>%  melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot)
    }   
    
    pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) + 
      scale_x_continuous(trans="log10", breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      scale_color_manual(name = NULL, values = group.colors) 
    
    if (res$Display_Ltot[1] == 1){ 
      pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Concentration of target–chaperone complex (in nM)")
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Concentration of target–chaperone complex (in nM)") 
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- pk + labs(x = "Total target concentration (in nM)", y = "Concentration of target–chaperone complex (in nM)") 
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Concentration of target–chaperone complex (in nM)") 
    }
    
    pk <- ggplotly(pk) # %>% layout(xaxis = list(exponentformat = "power"))
    print(pk)
  })
  
  
  # PLE and PL
  
  # Only ggplot (no plotly)
  # output$plot_PLE_PL <- renderPlot({
  #   res = res()
  #   
  #   if (res$Display_Ltot[1] == 1) { 
  #     pk <- res %>% select(Ltot, PLE, PL) %>% melt(id.vars = 'Ltot')  %>% dplyr::rename(xval = Ltot)
  #   } else if (res$Display_Ltot[1] == 0) { 
  #     pk <- res %>% select(Lfree, PLE, PL) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
  #   } else if (res$Display_Ltot[1] == 2) { 
  #     pk <- res %>% select(Ptot, PLE, PL) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot)
  #   } else if (res$Display_Ltot[1] == 3) { 
  #     pk <- res %>% select(Etot, PLE, PL) %>%  melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot)
  #   }   
  #   
  #   pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) +
  #       annotation_logticks(sides = "b") + scale_x_continuous(trans="log10", #limits=c(1e-1,NA),
  #                                                           breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                                                           labels = scales::trans_format("log10", scales::math_format(10^.x))
  #       ) + theme(text = element_text(size=14)) + scale_color_manual(name = NULL, values = group.colors)
  #   
  #   if (res$Display_Ltot[1] == 1){ 
  #     pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Concentration of ternary and\ntarget–ligand complex (in nM)")
  #   } else if (res$Display_Ltot[1] == 0) { 
  #     pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Concentration of ternary and\ntarget–ligand complex (in nM)") 
  #   } else if (res$Display_Ltot[1] == 2) { 
  #     pk <- pk + labs(x = "Total target concentration (in nM)", y = "Concentration of ternary and\ntarget–ligand complex (in nM)") 
  #   } else if (res$Display_Ltot[1] == 3) { 
  #     pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Concentration of ternary and\ntarget–ligand complex (in nM)") 
  #   }
  #   
  #   print(pk)
  # })
  
  # With plotly
  output$plot_PLE_PL <- renderPlotly({
    res = res()
    
    if (res$Display_Ltot[1] == 1) { 
      pk <- res %>% select(Ltot, PLE, PL) %>% dplyr::rename(TLC = PLE, TL = PL) %>% melt(id.vars = 'Ltot')  %>% dplyr::rename(xval = Ltot)
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- res %>% select(Lfree, PLE, PL) %>% dplyr::rename(TLC = PLE, TL = PL) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- res %>% select(Ptot, PLE, PL) %>% dplyr::rename(TLC = PLE, TL = PL) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot)
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- res %>% select(Etot, PLE, PL) %>% dplyr::rename(TLC = PLE, TL = PL) %>% melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot)
    }   
    
    pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) + 
      scale_x_continuous(trans="log10", breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      scale_color_manual(name = NULL, values = group.colors) 
    
    if (res$Display_Ltot[1] == 1){ 
      pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Concentration of ternary and\ntarget–ligand complex (in nM)")
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Concentration of ternary and\ntarget–ligand complex (in nM)") 
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- pk + labs(x = "Total target concentration (in nM)", y = "Concentration of ternary and\ntarget–ligand complex (in nM)") 
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Concentration of ternary and\ntarget–ligand complex (in nM)") 
    }
    
    pk <- ggplotly(pk) # %>% layout(xaxis = list(exponentformat = "power"))
    print(pk)
  })
  
  
  # PLE and EL
  
  # Only ggplot (no plotly)
  # output$plot_PLE_EL <- renderPlot({
  #   res = res()
  #   
  #   if (res$Display_Ltot[1] == 1) { 
  #     pk <- res %>% select(Ltot, PLE, EL) %>% melt(id.vars = 'Ltot')  %>% dplyr::rename(xval = Ltot)
  #   } else if (res$Display_Ltot[1] == 0) { 
  #     pk <- res %>% select(Lfree, PLE, EL) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
  #   } else if (res$Display_Ltot[1] == 2) { 
  #     pk <- res %>% select(Ptot, PLE, EL) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot)
  #   } else if (res$Display_Ltot[1] == 3) { 
  #     pk <- res %>% select(Etot, PLE, EL) %>%  melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot)
  #   }  
  #   
  #   pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) +
  #       annotation_logticks(sides = "b") + scale_x_continuous(trans="log10", #limits=c(1e-1,NA),
  #                                                           breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                                                           labels = scales::trans_format("log10", scales::math_format(10^.x))
  #       ) + theme(text = element_text(size=14)) + scale_color_manual(name = NULL, values = group.colors)#, labels = c("Ternary complex PLE", "Chaperone–Ligand EL")) 
  #     
  #   if (res$Display_Ltot[1] == 1){ 
  #     pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Concentration of ternary and\nchaperone–ligand complex (in nM)")
  #   } else if (res$Display_Ltot[1] == 0) { 
  #     pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Concentration of ternary and\nchaperone–ligand complex (in nM)") 
  #   } else if (res$Display_Ltot[1] == 2) { 
  #     pk <- pk + labs(x = "Total target concentration (in nM)", y = "Concentration of ternary and\nchaperone–ligand complex (in nM)") 
  #   } else if (res$Display_Ltot[1] == 3) { 
  #     pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Concentration of ternary and\nchaperone–ligand complex (in nM)") 
  #   }
  #   
  #   print(pk)
  # })
  # 
  
  # With plotly
  output$plot_PLE_EL <- renderPlotly({
    res = res()
    
    if (res$Display_Ltot[1] == 1) { 
      pk <- res %>% select(Ltot, PLE, EL) %>% dplyr::rename(TLC = PLE, CL = EL) %>% melt(id.vars = 'Ltot')  %>% dplyr::rename(xval = Ltot)
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- res %>% select(Lfree, PLE, EL) %>% dplyr::rename(TLC = PLE, CL = EL) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- res %>% select(Ptot, PLE, EL) %>% dplyr::rename(TLC = PLE, CL = EL) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot)
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- res %>% select(Etot, PLE, EL) %>% dplyr::rename(TLC = PLE, CL = EL) %>%  melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot)
    }  
    
    pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) + 
      scale_x_continuous(trans="log10", breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      scale_color_manual(name = NULL, values = group.colors) 
    
    if (res$Display_Ltot[1] == 1){ 
      pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Concentration of ternary and\nchaperone–ligand complex (in nM)")
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Concentration of ternary and\nchaperone–ligand complex (in nM)") 
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- pk + labs(x = "Total target concentration (in nM)", y = "Concentration of ternary and\nchaperone–ligand complex (in nM)") 
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Concentration of ternary and\nchaperone–ligand complex (in nM)") 
    }
    
    pk <- ggplotly(pk) # %>% layout(xaxis = list(exponentformat = "power"))
    print(pk)
  })
  
  
  # PLE and PE
  # With plotly
  output$plot_PLE_PE <- renderPlotly({
    res = res()
    
    if (res$Display_Ltot[1] == 1) { 
      pk <- res %>% select(Ltot, PLE, PE) %>% dplyr::rename(TLC = PLE, TC = PE) %>% melt(id.vars = 'Ltot')  %>% dplyr::rename(xval = Ltot)
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- res %>% select(Lfree, PLE, PE) %>% dplyr::rename(TLC = PLE, TC = PE) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- res %>% select(Ptot, PLE, PE) %>% dplyr::rename(TLC = PLE, TC = PE) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot)
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- res %>% select(Etot, PLE, PE) %>% dplyr::rename(TLC = PLE, TC = PE) %>%  melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot)
    }  
    
    pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) + 
      scale_x_continuous(trans="log10", breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      scale_color_manual(name = NULL, values = group.colors) 
    
    if (res$Display_Ltot[1] == 1){ 
      pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Concentration of ternary and\ntarget–chaperone complex (in nM)")
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Concentration of ternary and\ntarget–chaperone complex (in nM)") 
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- pk + labs(x = "Total target concentration (in nM)", y = "Concentration of ternary and\ntarget–chaperone complex (in nM)") 
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Concentration of ternary and\ntarget–chaperone complex (in nM)") 
    }
    
    pk <- ggplotly(pk) # %>% layout(xaxis = list(exponentformat = "power"))
    print(pk)
  })
  
  
  # Total target
  
  # Only ggplot (no plotly)
  # output$plot_tar_tot <- renderPlot({
  #   res = res()
  #   
  #   if (res$Display_Ltot[1] == 1) { 
  #     pk <- res %>% select(Ltot, Ptot) %>% melt(id.vars = 'Ltot')  %>% dplyr::rename(xval = Ltot)
  #   } else if (res$Display_Ltot[1] == 0) { 
  #     pk <- res %>% select(Lfree, Ptot) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
  #   } else if (res$Display_Ltot[1] == 2) { 
  #     pk <- res %>% select(Ptot) %>% mutate(Ptot_ = Ptot) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot) %>% mutate(variable = "Ptot")
  #   } else if (res$Display_Ltot[1] == 3) { 
  #     pk <- res %>% select(Etot, Ptot) %>%  melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot)
  #   }  
  #   
  #   pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) + 
  #     scale_x_continuous(trans="log10", breaks = scales::trans_breaks("log10", function(x) 10^x)) +
  #     scale_color_manual(name = NULL, values = group.colors) 
  #   
  #   if (res$Display_Ltot[1] == 1){ 
  #     pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Concentration of total target (in nM)")
  #   } else if (res$Display_Ltot[1] == 0) { 
  #     pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Concentration of total target (in nM)") 
  #   } else if (res$Display_Ltot[1] == 2) { 
  #     pk <- pk + labs(x = "Total target concentration (in nM)", y = "Concentration of total target (in nM)") 
  #   } else if (res$Display_Ltot[1] == 3) { 
  #     pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Concentration of total target (in nM)") 
  #   }
  #   
  #   print(pk)
  # })
  
  # With plotly
  output$plot_tar_tot <- renderPlotly({
    res = res()
    
    if (res$Display_Ltot[1] == 1) { 
      pk <- res %>% select(Ltot, Ptot) %>% dplyr::rename(Ttot = Ptot) %>% melt(id.vars = 'Ltot')  %>% dplyr::rename(xval = Ltot)
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- res %>% select(Lfree, Ptot) %>% dplyr::rename(Ttot = Ptot) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- res %>% select(Ptot) %>% mutate(Ttot = Ptot) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot) 
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- res %>% select(Etot, Ptot) %>% dplyr::rename(Ttot = Ptot) %>% melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot)
    }  
    
    pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) + 
     scale_x_continuous(trans="log10", breaks = scales::trans_breaks("log10", function(x) 10^x)) +
     scale_color_manual(name = NULL, values = group.colors) 
   
    if (res$Display_Ltot[1] == 1){ 
      pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Concentration of total target (in nM)")
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Concentration of total target (in nM)") 
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- pk + labs(x = "Total target concentration (in nM)", y = "Concentration of total target (in nM)") 
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Concentration of total target (in nM)") 
    }
    
    pk <- ggplotly(pk) # %>% layout(xaxis = list(exponentformat = "power"))
    print(pk)
  })
  
  
  # Free target
  
  # Only ggplot (no plotly)
  # output$plot_tar_free <- renderPlot({
  #   res = res()
  #   
  #   if (res$Display_Ltot[1] == 1) { 
  #     pk <- res %>% select(Ltot, P) %>% melt(id.vars = 'Ltot')  %>% dplyr::rename(xval = Ltot)
  #   } else if (res$Display_Ltot[1] == 0) { 
  #     pk <- res %>% select(Lfree, P) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
  #   } else if (res$Display_Ltot[1] == 2) { 
  #     pk <- res %>% select(Ptot, P) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot) 
  #   } else if (res$Display_Ltot[1] == 3) { 
  #     pk <- res %>% select(Etot, P) %>%  melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot)
  #   }  
  #  
  #   pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) +
  #       annotation_logticks(sides = "b") + scale_x_continuous(trans="log10", #limits=c(1e-1,NA),
  #                                                           breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                                                           labels = scales::trans_format("log10", scales::math_format(10^.x))
  #       ) + theme(text = element_text(size=14)) + scale_color_manual(name = NULL, values = group.colors, labels = c("Free target"))
  #  
  #   if (res$Display_Ltot[1] == 1){ 
  #     pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Concentration of free target (in nM)")
  #   } else if (res$Display_Ltot[1] == 0) { 
  #     pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Concentration of free target (in nM)") 
  #   } else if (res$Display_Ltot[1] == 2) { 
  #     pk <- pk + labs(x = "Total target concentration (in nM)", y = "Concentration of free target (in nM)") 
  #   } else if (res$Display_Ltot[1] == 3) { 
  #     pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Concentration of free target (in nM)") 
  #   }
  #   
  #   print(pk)
  # })
  
  # With plotly
  output$plot_tar_free <- renderPlotly({
    res = res()
    
    if (res$Display_Ltot[1] == 1) { 
      pk <- res %>% select(Ltot, Pfree) %>% dplyr::rename(Tfree = Pfree) %>% melt(id.vars = 'Ltot')  %>% dplyr::rename(xval = Ltot)
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- res %>% select(Lfree, Pfree) %>% dplyr::rename(Tfree = Pfree) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- res %>% select(Ptot, Pfree) %>% dplyr::rename(Tfree = Pfree) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot) 
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- res %>% select(Etot, Pfree) %>% dplyr::rename(Tfree = Pfree) %>% melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot)
    }  
    
    pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) + 
      scale_x_continuous(trans="log10", breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      scale_color_manual(name = NULL, values = group.colors) 
    
    if (res$Display_Ltot[1] == 1){ 
      pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Concentration of free target (in nM)")
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Concentration of free target (in nM)") 
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- pk + labs(x = "Total target concentration (in nM)", y = "Concentration of free target (in nM)") 
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Concentration of free target (in nM)") 
    }
    
    pk <- ggplotly(pk) # %>% layout(xaxis = list(exponentformat = "power"))
    print(pk)
  })
 
  
  # TO / Percent target bound

  # Only ggplot (no plotly)
  # output$plot_TO <- renderPlot({
  #   res = res()
  #   
  #   if (res$Display_Ltot[1] == 1) { 
  #     pk <- res %>% select(Ltot, P_bnd) %>% melt(id.vars = 'Ltot')  %>% dplyr::rename(xval = Ltot)
  #   } else if (res$Display_Ltot[1] == 0) { 
  #     pk <- res %>% select(Lfree, P_bnd) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
  #   } else if (res$Display_Ltot[1] == 2) { 
  #     pk <- res %>% select(Ptot, P_bnd) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot) 
  #   } else if (res$Display_Ltot[1] == 3) { 
  #     pk <- res %>% select(Etot, P_bnd) %>%  melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot)
  #   }  
  #   
  #   pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) +
  #   annotation_logticks(sides = "b") + scale_x_continuous(trans="log10", #limits=c(1e-1,NA),
  #                                                           breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                                                           labels = scales::trans_format("log10", scales::math_format(10^.x))
  #   ) + theme(text = element_text(size=14)) + scale_color_manual(name = NULL, values = group.colors, labels = c("Percent target bound"))
  #     
  #   if (res$Display_Ltot[1] == 1){ 
  #     pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Percent of target bound / target occupacy (in Percent)")
  #   } else if (res$Display_Ltot[1] == 0) { 
  #     pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Percent of target bound / target occupacy (in Percent)") 
  #   } else if (res$Display_Ltot[1] == 2) { 
  #     pk <- pk + labs(x = "Total target concentration (in nM)", y = "Percent of target bound / target occupacy (in Percent)") 
  #   } else if (res$Display_Ltot[1] == 3) { 
  #     pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Percent of target bound / target occupacy (in Percent)") 
  #   }
  #   
  #   print(pk)
  # })
  
  # With plotly
  output$plot_TO <- renderPlotly({
    res = res()
    Ltot_EC50_Tbound <- as.numeric(res$Ltot_EC50_Tbound)
    
    if (res$Display_Ltot[1] == 1) { 
      pk <- res %>% select(Ltot, P_bnd) %>% dplyr::rename(T_bound = P_bnd) %>% melt(id.vars = 'Ltot')  %>% dplyr::rename(xval = Ltot)
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- res %>% select(Lfree, P_bnd) %>% dplyr::rename(T_bound = P_bnd) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- res %>% select(Ptot, P_bnd) %>% dplyr::rename(T_bound = P_bnd) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot) 
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- res %>% select(Etot, P_bnd) %>% dplyr::rename(T_bound = P_bnd) %>% melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot)
    }  
    
    pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) + 
      scale_x_continuous(trans="log10", breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      scale_color_manual(name = NULL, values = group.colors) 
    
    if (res$Display_Ltot[1] == 1){ 
      pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Percent of target bound / target occupacy (in Percent)") +
        geom_vline(xintercept = Ltot_EC50_Tbound, linetype="dashed", color = "black", size=1) + annotate("text", x = Ltot_EC50_Tbound, y = res$P_bnd %>% min(), hjust = 0, vjust = -0.25, angle = 90, label = "EC50", size = 4, parse = TRUE)
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Percent of target bound / target occupacy (in Percent)")         
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- pk + labs(x = "Total target concentration (in nM)", y = "Percent of target bound / target occupacy (in Percent)") 
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Percent of target bound / target occupacy (in Percent)") 
    }
    
    pk <- ggplotly(pk) # %>% layout(xaxis = list(exponentformat = "power"))
    print(pk)
  })
  
  
  # Percent target & chaperone bound
  
  # With plotly
  output$plot_Tbnd_Cbnd <- renderPlotly({
    res = res()
    
    if (res$Display_Ltot[1] == 1) { 
      pk <- res %>% select(Ltot, P_bnd, E_bnd) %>% dplyr::rename(T_bound = P_bnd, C_bound = E_bnd) %>% melt(id.vars = 'Ltot')  %>% dplyr::rename(xval = Ltot)
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- res %>% select(Lfree, P_bnd, E_bnd) %>% dplyr::rename(T_bound = P_bnd, C_bound = E_bnd) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- res %>% select(Ptot, P_bnd, E_bnd) %>% dplyr::rename(T_bound = P_bnd, C_bound = E_bnd) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot) 
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- res %>% select(Etot, P_bnd, E_bnd) %>% dplyr::rename(T_bound = P_bnd, C_bound = E_bnd) %>% melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot)
    }  
    
    pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) + 
      scale_x_continuous(trans="log10", breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      scale_color_manual(name = NULL, values = group.colors) 
    
    if (res$Display_Ltot[1] == 1){ 
      pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "% of target and chaperone bound")
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "% of target and chaperone bound") 
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- pk + labs(x = "Total target concentration (in nM)", y = "% of target and chaperone bound") 
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "% of target and chaperone bound") 
    }
    
    pk <- ggplotly(pk) # %>% layout(xaxis = list(exponentformat = "power"))
    print(pk)
  })
  
  # Target bound in ternary complex
  
  # Only ggplot (no plotly)
  #   output$plot_PLE_P <- renderPlot({
  #     res = res()
  #     
  #     if (res$Display_Ltot[1] == 1) { 
  #       pk <- res %>% select(Ltot, PLE_Ptot) %>% melt(id.vars = 'Ltot')  %>% dplyr::rename(xval = Ltot)
  #     } else if (res$Display_Ltot[1] == 0) { 
  #       pk <- res %>% select(Lfree, PLE_Ptot) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
  #     } else if (res$Display_Ltot[1] == 2) { 
  #       pk <- res %>% select(Ptot, PLE_Ptot) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot) 
  #     } else if (res$Display_Ltot[1] == 3) { 
  #       pk <- res %>% select(Etot, PLE_Ptot) %>%  melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot)
  #     }  
  #     
  #     pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) +
  #     annotation_logticks(sides = "b") + scale_x_continuous(trans="log10", #limits=c(1e-1,NA),
  #                                                             breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                                                             labels = scales::trans_format("log10", scales::math_format(10^.x))
  #     ) + theme(text = element_text(size=14)) + scale_color_manual(name = NULL, values = group.colors, labels = c("Percent of PLE/Ptot")) 
  #     
  #     if (res$Display_Ltot[1] == 1){ 
  #       pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Percent of ternary complex in total target")
  #     } else if (res$Display_Ltot[1] == 0) { 
  #       pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Percent of ternary complex in total target") 
  #     } else if (res$Display_Ltot[1] == 2) { 
  #       pk <- pk + labs(x = "Total target concentration (in nM)", y = "Percent of ternary complex in total target") 
  #     } else if (res$Display_Ltot[1] == 3) { 
  #       pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Percent of ternary complex in total target") 
  #     }
  #       
  #     print(pk)
  #   })
  
  # With plotly
  output$plot_PLE_P <- renderPlotly({
    res = res()
    
    if (res$Display_Ltot[1] == 1) { 
      pk <- res %>% select(Ltot, PLE_Ptot) %>% dplyr::rename(TLC_Ttot = PLE_Ptot) %>% melt(id.vars = 'Ltot')  %>% dplyr::rename(xval = Ltot)
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- res %>% select(Lfree, PLE_Ptot) %>% dplyr::rename(TLC_Ttot = PLE_Ptot) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- res %>% select(Ptot, PLE_Ptot) %>% dplyr::rename(TLC_Ttot = PLE_Ptot) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot) 
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- res %>% select(Etot, PLE_Ptot) %>% dplyr::rename(TLC_Ttot = PLE_Ptot) %>% melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot)
    }  
    
    pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) + 
      scale_x_continuous(trans="log10", breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      scale_color_manual(name = NULL, values = group.colors) 
    
    if (res$Display_Ltot[1] == 1){ 
      pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Percent of ternary complex in total target")
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Percent of ternary complex in total target") 
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- pk + labs(x = "Total target concentration (in nM)", y = "Percent of ternary complex in total target") 
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Percent of ternary complex in total target") 
    }
    
    pk <- ggplotly(pk) # %>% layout(xaxis = list(exponentformat = "power"))
    print(pk)
  })
  
  
  # Total Ligand
  
  # Only ggplot (no plotly)
  # output$plot_lig_tot <- renderPlot({
  #   res = res()
  #   TF50 <- res$TF50[1]
  #   TI50 <- res$TI50[1]
  #   Ltotmax <- res$Ltotmax[1]
  #   
  #   pk <- res %>% select(Lfree, Ltot) %>% 
  #     melt(id.vars = 'Lfree') %>% 
  #     ggplot(aes(x =  Lfree, y= value, color = variable)) + geom_line(alpha=1, size=1.25) +
  #     annotation_logticks(sides = "b") + scale_x_continuous(trans="log10", #limits=c(1e-1,NA),
  #                                                           breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                                                           labels = scales::trans_format("log10", scales::math_format(10^.x))
  #     ) +
  #     annotation_logticks(sides = "l") + scale_y_continuous(trans="log10", #limits=c(1e-1,NA),
  #                                                           breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                                                           labels = scales::trans_format("log10", scales::math_format(10^.x))
  #     ) +
  #     geom_hline(yintercept = as.numeric(TF50), linetype = "dashed", size = .75) + annotate("text", y = as.numeric(TF50), x = res$Lfree %>% max(), vjust = -0.15, hjust = 0.75, label = "TF[50]", size = 6, parse = TRUE) + 
  #     geom_hline(yintercept = as.numeric(TI50), linetype = "dashed", size = .75) + annotate("text", y = as.numeric(TI50), x = res$Lfree %>% max(), vjust = -0.15, hjust = 0.75, label = "TI[50]", size = 6, parse = TRUE) + 
  #     geom_hline(yintercept = as.numeric(Ltotmax), linetype = "dashed", size = .75) + annotate("text", y = as.numeric(Ltotmax), x = res$Lfree %>% max(), vjust = -0.15, hjust = 0.75, label = "L['tot,max']", size = 6, parse = TRUE) + 
  #     labs(x = "Free ligand concentration (estimated)  (in nM)", y = "Total ligand concentration (in nM)") + theme(text = element_text(size=14)) +
  #     scale_color_manual(name = NULL, values = group.colors, labels = c("Total ligand"))
  #   
  #   print(pk)
  # })

  # With plotly
  output$plot_lig_tot <- renderPlotly({
    res = res()
    TF50 <- res$TF50[1]
    TI50 <- res$TI50[1]
    Ltotmax <- res$Ltotmax[1]
    
    pk <- res %>% select(Lfree, Ltot) %>% 
      melt(id.vars = 'Lfree') %>% 
      ggplot(aes(x =  Lfree, y = value, color = variable)) + geom_line(alpha=1, size=1.25) +
      scale_x_continuous(trans="log10", breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      scale_y_continuous(trans="log10", breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      scale_color_manual(name = NULL, values = group.colors) +
      geom_hline(yintercept = as.numeric(TF50), linetype = "dashed", size = .75) + annotate("text", y = as.numeric(TF50)*0.6, x = res$Lfree %>% max()*0.7, vjust = -0.15, hjust = 0.75, label = "TF50", size = 4, parse = TRUE) + 
      geom_hline(yintercept = as.numeric(TI50), linetype = "dashed", size = .75) + annotate("text", y = as.numeric(TI50)*0.6, x = res$Lfree %>% max()*0.7, vjust = -0.15, hjust = 0.75, label = "TI50", size = 4, parse = TRUE) + 
      geom_hline(yintercept = as.numeric(Ltotmax), linetype = "dashed", size = .75) + annotate("text", y = as.numeric(Ltotmax)*0.6, x = res$Lfree %>% max()*0.5, vjust = -0.15, hjust = 0.75, label = "Ltot,max", size = 4, parse = TRUE) + 
      labs(x = "Free ligand concentration (estimated)  (in nM)", y = "Total ligand concentration (in nM)") 
    
    pk <- ggplotly(pk) # %>% layout(xaxis = list(exponentformat = "power"))
    print(pk)
  })
  
    
  # Percent ligand bound
  
  # Only ggplot (no plotly)
  # output$plot_lig_bound <- renderPlot({
  #   
  #   res = res()
  #   if (res$Display_Ltot[1] == 1) { 
  #     pk <- res %>% select(Ltot, L_bnd) %>% melt(id.vars = 'Ltot') %>% dplyr::rename(xval = Ltot)
  #   } else if (res$Display_Ltot[1] == 0) { 
  #     pk <- res %>% select(Lfree, L_bnd) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
  #   } else if (res$Display_Ltot[1] == 2) { 
  #     pk <- res %>% select(Ptot, L_bnd) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot)
  #   } else if (res$Display_Ltot[1] == 3) { 
  #     pk <- res %>% select(Etot, L_bnd) %>% melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot)
  #   }         
  #   
  #   pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable))  + geom_line(alpha=1, size=1.25) +
  #       annotation_logticks(sides = "b") + scale_x_continuous(trans="log10", #limits=c(1e-1,NA),
  #                                                             breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                                                             labels = scales::trans_format("log10", scales::math_format(10^.x))
  #       ) +
  #       scale_color_manual(name = NULL, values = group.colors, labels = c("Percent ligand bound")) + theme(text = element_text(size=14)) 
  #   
  #   if (res$Display_Ltot[1] == 1){ 
  #     pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Ligand bound (in Percent)") 
  #   } else if (res$Display_Ltot[1] == 0) { 
  #     pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Ligand bound (in Percent)") 
  #   } else if (res$Display_Ltot[1] == 2) { 
  #     pk <- pk + labs(x = "Total target concentration (in nM)", y = "Ligand bound (in Percent)") 
  #   } else if (res$Display_Ltot[1] == 3) { 
  #     pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Ligand bound (in Percent)") 
  #   }
  # 
  #   print(pk)
  # })
  
  # With plotly
  output$plot_lig_bound <- renderPlotly({
    
    res = res()
    if (res$Display_Ltot[1] == 1) { 
      pk <- res %>% select(Ltot, L_bnd) %>% dplyr::rename(L_bound = L_bnd) %>% melt(id.vars = 'Ltot') %>% dplyr::rename(xval = Ltot)
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- res %>% select(Lfree, L_bnd) %>% dplyr::rename(L_bound = L_bnd) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- res %>% select(Ptot, L_bnd) %>% dplyr::rename(L_bound = L_bnd) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot)
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- res %>% select(Etot, L_bnd) %>% dplyr::rename(L_bound = L_bnd) %>% melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot)
    }         
    
    pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) + 
      scale_x_continuous(trans="log10", breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      scale_color_manual(name = NULL, values = group.colors) 
    
    if (res$Display_Ltot[1] == 1){ 
      pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Ligand bound (in Percent)") 
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Ligand bound (in Percent)") 
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- pk + labs(x = "Total target concentration (in nM)", y = "Ligand bound (in Percent)") 
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Ligand bound (in Percent)") 
    }
    
    pk <- ggplotly(pk) # %>% layout(xaxis = list(exponentformat = "power"))
    print(pk)
  })
  
    
  # Total chaperone
  
  # Only ggplot (no plotly)
  # output$plot_chap_tot <- renderPlot({
  #   res = res()
  #   
  #   if (res$Display_Ltot[1] == 1) { 
  #     pk <- res %>% select(Ltot, Etot) %>% melt(id.vars = 'Ltot') %>% dplyr::rename(xval = Ltot)
  #   } else if (res$Display_Ltot[1] == 0) { 
  #     pk <- res %>% select(Lfree, Etot) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
  #   } else if (res$Display_Ltot[1] == 2) { 
  #     pk <- res %>% select(Ptot, Etot) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot)
  #   } else if (res$Display_Ltot[1] == 3) { 
  #     pk <- res %>% select(Etot) %>% mutate(Etot_ = Etot) %>% melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot) %>% mutate(variable = "Etot")
  #   }         
  #   
  #   pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) + 
  #     annotation_logticks(sides = "b") + scale_x_continuous(trans="log10", #limits=c(1e-1,NA),
  #                                                           breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                                                           labels = scales::trans_format("log10", scales::math_format(10^.x))
  #       ) + theme(text = element_text(size=14)) + scale_color_manual(name = NULL, values = group.colors, labels = c("Total chaperone")) 
  #   
  #   if (res$Display_Ltot[1] == 1){ 
  #     pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Concentration of total chaperone (in nM)")
  #   } else if (res$Display_Ltot[1] == 0) { 
  #     pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Concentration of total chaperone (in nM)") 
  #   } else if (res$Display_Ltot[1] == 2) { 
  #     pk <- pk + labs(x = "Total target concentration (in nM)", y = "Concentration of total chaperone (in nM)") 
  #   } else if (res$Display_Ltot[1] == 3) { 
  #     pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Concentration of total chaperone (in nM)") 
  #   }
  #   
  #   print(pk)
  # })
  
  # With plotly
  output$plot_chap_tot <- renderPlotly({
    res = res()
    
    if (res$Display_Ltot[1] == 1) { 
      pk <- res %>% select(Ltot, Etot) %>% dplyr::rename(Ctot = Etot) %>% melt(id.vars = 'Ltot') %>% dplyr::rename(xval = Ltot)
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- res %>% select(Lfree, Etot) %>% dplyr::rename(Ctot = Etot) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- res %>% select(Ptot, Etot) %>% dplyr::rename(Ctot = Etot) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot)
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- res %>% select(Etot) %>% mutate(Ctot = Etot) %>% melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot)
    }         
    
    pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) + 
      scale_x_continuous(trans="log10", breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      scale_color_manual(name = NULL, values = group.colors) 
    
    if (res$Display_Ltot[1] == 1){ 
      pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Concentration of total chaperone (in nM)")
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Concentration of total chaperone (in nM)") 
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- pk + labs(x = "Total target concentration (in nM)", y = "Concentration of total chaperone (in nM)") 
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Concentration of total chaperone (in nM)") 
    }
    
    pk <- ggplotly(pk) # %>% layout(xaxis = list(exponentformat = "power"))
    print(pk)
  })
  
  
  # Free chaperone
  
  # Only ggplot (no plotly)
  # output$plot_chap_free <- renderPlot({
  #   res = res()
  #   
  #   if (res$Display_Ltot[1] == 1) { 
  #     pk <- res %>% select(Ltot, E) %>% melt(id.vars = 'Ltot') %>% dplyr::rename(xval = Ltot)
  #   } else if (res$Display_Ltot[1] == 0) { 
  #     pk <- res %>% select(Lfree, E) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
  #   } else if (res$Display_Ltot[1] == 2) { 
  #     pk <- res %>% select(Ptot, E) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot)
  #   } else if (res$Display_Ltot[1] == 3) { 
  #     pk <- res %>% select(Etot, E) %>% melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot) 
  #   }         
  #   
  #   pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) + 
  #   annotation_logticks(sides = "b") + scale_x_continuous(trans="log10", #limits=c(1e-1,NA),
  #                                                             breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                                                             labels = scales::trans_format("log10", scales::math_format(10^.x))
  #   ) + theme(text = element_text(size=14)) + scale_color_manual(name = NULL, values = group.colors, labels = c("Free chaperone")) 
  #   
  #   if (res$Display_Ltot[1] == 1){ 
  #     pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Concentration of free chaperone (in nM)")
  #   } else if (res$Display_Ltot[1] == 0) { 
  #     pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Concentration of free chaperone (in nM)") 
  #   } else if (res$Display_Ltot[1] == 2) { 
  #     pk <- pk + labs(x = "Total target concentration (in nM)", y = "Concentration of free chaperone (in nM)") 
  #   } else if (res$Display_Ltot[1] == 3) { 
  #     pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Concentration of free chaperone (in nM)") 
  #   }
  #   
  #   print(pk)
  # })
  
  # With plotly
  output$plot_chap_free <- renderPlotly({
    res = res()
    
    if (res$Display_Ltot[1] == 1) { 
      pk <- res %>% select(Ltot, Efree) %>% dplyr::rename(Cfree = Efree) %>% melt(id.vars = 'Ltot') %>% dplyr::rename(xval = Ltot)
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- res %>% select(Lfree, Efree) %>% dplyr::rename(Cfree = Efree) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- res %>% select(Ptot, Efree) %>% dplyr::rename(Cfree = Efree) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot)
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- res %>% select(Etot, Efree) %>% dplyr::rename(Cfree = Efree) %>% melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot) 
    }         
    
    pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) + 
      scale_x_continuous(trans="log10", breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      scale_color_manual(name = NULL, values = group.colors) 
    
    if (res$Display_Ltot[1] == 1){ 
      pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Concentration of free chaperone (in nM)")
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Concentration of free chaperone (in nM)") 
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- pk + labs(x = "Total target concentration (in nM)", y = "Concentration of free chaperone (in nM)") 
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Concentration of free chaperone (in nM)") 
    }
    
    pk <- ggplotly(pk) # %>% layout(xaxis = list(exponentformat = "power"))
    print(pk)
  })
  
  
  # Percent chaperone bound
  
  # Only ggplot (no plotly)
  # output$plot_chap_bound <- renderPlot({
  #   res = res()
  #   
  #   if (res$Display_Ltot[1] == 1) { 
  #     pk <- res %>% select(Ltot, E_bnd) %>% melt(id.vars = 'Ltot') %>% dplyr::rename(xval = Ltot)
  #   } else if (res$Display_Ltot[1] == 0) { 
  #     pk <- res %>% select(Lfree, E_bnd) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
  #   } else if (res$Display_Ltot[1] == 2) { 
  #     pk <- res %>% select(Ptot, E_bnd) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot)
  #   } else if (res$Display_Ltot[1] == 3) { 
  #     pk <- res %>% select(Etot, E_bnd) %>% melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot) 
  #   }         
  #   
  #   pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) + 
  #     annotation_logticks(sides = "b") + scale_x_continuous(trans="log10", #limits=c(1e-1,NA),
  #                                                           breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                                                           labels = scales::trans_format("log10", scales::math_format(10^.x))
  #     ) + theme(text = element_text(size=14)) + scale_color_manual(name = NULL, values = group.colors, labels = c("Percent chaperone bound")) 
  #   
  #   if (res$Display_Ltot[1] == 1){ 
  #     pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Chaperone bound (in Percent)")
  #   } else if (res$Display_Ltot[1] == 0) { 
  #     pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Chaperone bound (in Percent)") 
  #   } else if (res$Display_Ltot[1] == 2) { 
  #     pk <- pk + labs(x = "Total target concentration (in nM)", y = "Chaperone bound (in Percent)") 
  #   } else if (res$Display_Ltot[1] == 3) { 
  #     pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Chaperone bound (in Percent)") 
  #   }
  #   
  #   print(pk)
  # })
  
  # With plotly
  output$plot_chap_bound <- renderPlotly({
    res = res()
    Ltot_EC50_Cbound <- as.numeric(res$Ltot_EC50_Cbound)
    
    if (res$Display_Ltot[1] == 1) { 
      pk <- res %>% select(Ltot, E_bnd) %>% dplyr::rename(C_bound = E_bnd) %>% melt(id.vars = 'Ltot') %>% dplyr::rename(xval = Ltot)
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- res %>% select(Lfree, E_bnd) %>% dplyr::rename(C_bound = E_bnd) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- res %>% select(Ptot, E_bnd) %>% dplyr::rename(C_bound = E_bnd) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot)
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- res %>% select(Etot, E_bnd) %>% dplyr::rename(C_bound = E_bnd) %>% melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot) 
    }         
    
    pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) + 
      scale_x_continuous(trans="log10", breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      scale_color_manual(name = NULL, values = group.colors) 
    
    if (res$Display_Ltot[1] == 1){ 
      pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Chaperone bound (in Percent)") +
        geom_vline(xintercept = Ltot_EC50_Cbound, linetype="dashed", color = "black", size=1) + annotate("text", x = Ltot_EC50_Cbound, y = res$P_bnd %>% min(), hjust = 0, vjust = -0.25, angle = 90, label = "EC50", size = 4, parse = TRUE)
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Chaperone bound (in Percent)") 
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- pk + labs(x = "Total target concentration (in nM)", y = "Chaperone bound (in Percent)") 
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Chaperone bound (in Percent)") 
    }
    
    pk <- ggplotly(pk) # %>% layout(xaxis = list(exponentformat = "power"))
    print(pk)
  })
  
  # Percent target & chaperone bound
  
  # With plotly
  output$plot_Cbnd_Tbnd <- renderPlotly({
    res = res()
    
    if (res$Display_Ltot[1] == 1) { 
      pk <- res %>% select(Ltot, P_bnd, E_bnd) %>% dplyr::rename(T_bound = P_bnd, C_bound = E_bnd) %>% melt(id.vars = 'Ltot')  %>% dplyr::rename(xval = Ltot)
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- res %>% select(Lfree, P_bnd, E_bnd) %>% dplyr::rename(T_bound = P_bnd, C_bound = E_bnd) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- res %>% select(Ptot, P_bnd, E_bnd) %>% dplyr::rename(T_bound = P_bnd, C_bound = E_bnd) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot) 
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- res %>% select(Etot, P_bnd, E_bnd) %>% dplyr::rename(T_bound = P_bnd, C_bound = E_bnd) %>% melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot)
    }  
    
    pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) + 
      scale_x_continuous(trans="log10", breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      scale_color_manual(name = NULL, values = group.colors) 
    
    if (res$Display_Ltot[1] == 1){ 
      pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "% of target and chaperone bound")
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "% of target and chaperone bound") 
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- pk + labs(x = "Total target concentration (in nM)", y = "% of target and chaperone bound") 
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "% of target and chaperone bound") 
    }
    
    pk <- ggplotly(pk) # %>% layout(xaxis = list(exponentformat = "power"))
    print(pk)
  })
  
  # Chaperone bound in ternary complex
  
  # Only ggplot (no plotly)
  # output$plot_PLE_E <- renderPlot({
  #   res = res()
  #   
  #   if (res$Display_Ltot[1] == 1) { 
  #     pk <- res %>% select(Ltot, PLE_Etot) %>% melt(id.vars = 'Ltot') %>% dplyr::rename(xval = Ltot)
  #   } else if (res$Display_Ltot[1] == 0) { 
  #     pk <- res %>% select(Lfree, PLE_Etot) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
  #   } else if (res$Display_Ltot[1] == 2) { 
  #     pk <- res %>% select(Ptot, PLE_Etot) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot)
  #   } else if (res$Display_Ltot[1] == 3) { 
  #     pk <- res %>% select(Etot, PLE_Etot) %>% melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot) 
  #   }         
  #   
  #   pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) + 
  #     annotation_logticks(sides = "b") + scale_x_continuous(trans="log10", #limits=c(1e-1,NA),
  #                                                           breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                                                           labels = scales::trans_format("log10", scales::math_format(10^.x))
  #     ) + theme(text = element_text(size=14)) + scale_color_manual(name = NULL, values = group.colors, labels = c("Percent of PLE/Etot")) 
  #   
  #   if (res$Display_Ltot[1] == 1){ 
  #     pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Percent of ternary complex in total chaperone")
  #   } else if (res$Display_Ltot[1] == 0) { 
  #     pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Percent of ternary complex in total chaperone") 
  #   } else if (res$Display_Ltot[1] == 2) { 
  #     pk <- pk + labs(x = "Total target concentration (in nM)", y = "Percent of ternary complex in total chaperone") 
  #   } else if (res$Display_Ltot[1] == 3) { 
  #     pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Percent of ternary complex in total chaperone") 
  #   }
  #   
  #   print(pk)
  # })
  
  # With plotly
  output$plot_PLE_E <- renderPlotly({
    res = res()
    
    if (res$Display_Ltot[1] == 1) { 
      pk <- res %>% select(Ltot, PLE_Etot) %>% dplyr::rename(TLC_Ctot = PLE_Etot) %>% melt(id.vars = 'Ltot') %>% dplyr::rename(xval = Ltot)
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- res %>% select(Lfree, PLE_Etot) %>% dplyr::rename(TLC_Ctot = PLE_Etot) %>% melt(id.vars = 'Lfree') %>% dplyr::rename(xval = Lfree)
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- res %>% select(Ptot, PLE_Etot) %>% dplyr::rename(TLC_Ctot = PLE_Etot) %>% melt(id.vars = 'Ptot') %>% dplyr::rename(xval = Ptot)
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- res %>% select(Etot, PLE_Etot) %>% dplyr::rename(TLC_Ctot = PLE_Etot) %>% melt(id.vars = 'Etot') %>% dplyr::rename(xval = Etot) 
    }         
    
    pk <- pk %>% ggplot(aes(x = xval, y = value, color = variable)) + geom_line(alpha=1, size=1.25) + 
      scale_x_continuous(trans="log10", breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      scale_color_manual(name = NULL, values = group.colors) 
    
    if (res$Display_Ltot[1] == 1){ 
      pk <- pk + labs(x = "Total ligand concentration (in nM)", y = "Percent of ternary complex in total chaperone")
    } else if (res$Display_Ltot[1] == 0) { 
      pk <- pk + labs(x = "Free ligand concentration (estimated) (in nM)", y = "Percent of ternary complex in total chaperone") 
    } else if (res$Display_Ltot[1] == 2) { 
      pk <- pk + labs(x = "Total target concentration (in nM)", y = "Percent of ternary complex in total chaperone") 
    } else if (res$Display_Ltot[1] == 3) { 
      pk <- pk + labs(x = "Total chaperone concentration (in nM)", y = "Percent of ternary complex in total chaperone") 
    }
    
    pk <- ggplotly(pk) # %>% layout(xaxis = list(exponentformat = "power"))
    print(pk)
  })
  
  
  # Tracer and Tracer bound
  
  # Only ggplot (no plotly)
  # output$plot_EC_E <- renderPlot({
  #   res = res()
  #   EC_Lsmall <- res$EC_Lsmall[1]
  #   
  #   if (res$Display_Ltot[1] == 1){
  #     pk <- res %>% select(Ltot, E, EC) %>% melt(id.vars = 'Ltot') %>% 
  #       ggplot(aes(x = Ltot, y = value, color = variable)) + geom_line(alpha=1, size=1.25) +
  #       annotation_logticks(sides = "b") + scale_x_continuous(trans="log10", #limits=c(1e-1,NA),
  #                                                             breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                                                             labels = scales::trans_format("log10", scales::math_format(10^.x))
  #       ) +
  #       geom_hline(yintercept = EC_Lsmall, linetype = "dashed", size = .75) + annotate("text", x = res$Ltot %>% last(), y = as.numeric(EC_Lsmall), hjust = 0.75, vjust = 1.25, label = "EC['lim']", size = 5, parse = TRUE) + 
  #       labs(x = "Total ligand concentration (in nM)", y = "Concentration of tracer and complexes (in nM)") + theme(text = element_text(size=14)) +
  #       scale_color_manual(name = NULL, values = group.colors, labels = c("E", "EC")) 
  #   } else { 
  #     pk <- res %>% select(Lfree, E, EC) %>% melt(id.vars = 'Lfree') %>% 
  #       ggplot(aes(x = Lfree, y = value, color = variable)) + geom_line(alpha=1, size=1.25) +
  #       annotation_logticks(sides = "b") + scale_x_continuous(trans="log10", #limits=c(1e-1,NA),
  #                                                             breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                                                             labels = scales::trans_format("log10", scales::math_format(10^.x))
  #       ) +
  #       geom_hline(yintercept = EC_Lsmall, linetype = "dashed", size = .75) + annotate("text", x = res$Lfree %>% last(), y = as.numeric(EC_Lsmall), hjust = 0.75, vjust = 1.25, label = "EC['lim']", size = 5, parse = TRUE) + 
  #       labs(x = "Free ligand concentration (estimated) (in nM)", y = "Concentration of tracer and complexes (in nM)") + theme(text = element_text(size=14)) +
  #       scale_color_manual(name = NULL, values = group.colors, labels = c("EC", "C")) 
  #   }
  #   
  #   print(pk)
  # })
  
  # plotly
  output$plot_EC_E <- renderPlotly({
    res = res()
    EC_Lsmall <- res$EC_Lsmall[1]
    
    if (res$Display_Ltot[1] == 1){
      pk <- res %>% select(Ltot, Cfree, EC) %>% dplyr::rename(Trfree = Cfree, CTr = EC) %>% melt(id.vars = 'Ltot') %>% 
        ggplot(aes(x = Ltot, y = value, color = variable)) + geom_line(alpha=1, size=1.25) +
        scale_x_continuous(trans="log10", breaks = scales::trans_breaks("log10", function(x) 10^x)) +
        geom_hline(yintercept = EC_Lsmall, linetype = "dashed", size = .75) + annotate("text", x = res$Ltot %>% last()*0.5, y = as.numeric(EC_Lsmall)*1.1, hjust = 0.75, vjust = 0.25, label = "C-Trlim", size = 4, parse = TRUE) + 
        labs(x = "Total ligand concentration (in nM)", y = "Concentration of tracer and complexes (in nM)") +
        scale_color_manual(name = NULL, values = group.colors) 
    } else { 
      pk <- res %>% select(Lfree, Cfree, EC) %>% dplyr::rename(Trfree = Cfree, CTr = EC) %>% melt(id.vars = 'Lfree') %>% 
        ggplot(aes(x = Lfree, y = value, color = variable)) + geom_line(alpha=1, size=1.25) +
        scale_x_continuous(trans="log10", breaks = scales::trans_breaks("log10", function(x) 10^x)) +
        geom_hline(yintercept = EC_Lsmall, linetype = "dashed", size = .75) + annotate("text", x = res$Ltot %>% last()*0.5, y = as.numeric(EC_Lsmall)*1.1, hjust = 0.75, vjust = 0.25, label = "C-Trlim", size = 4, parse = TRUE) + 
        labs(x = "Free ligand concentration (estimated) (in nM)", y = "Concentration of tracer and complexes (in nM)") +
        scale_color_manual(name = NULL, values = group.colors) 
    }
    
    pk <- ggplotly(pk) # %>% layout(xaxis = list(exponentformat = "power"))
    print(pk)
  })
  
  
  # Tracer, tracer bound and complexes
  
  # Only ggplot (no plotly)
  # output$plot_EC_E_E_Ctot <- renderPlot({
  #   res = res()
  #   EC_Lsmall <- res$EC_Lsmall[1]
  #   
  #   if (res$Display_Ltot[1] == 1){
  #     pk <- res %>% select(Ltot, E, EC, Cfree, Ctot) %>% melt(id.vars = 'Ltot') %>% 
  #       ggplot(aes(x = Ltot, y = value, color = variable)) + geom_line(alpha=1, size=1.25) +
  #       annotation_logticks(sides = "b") + scale_x_continuous(trans="log10", #limits=c(1e-1,NA),
  #                                                             breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                                                             labels = scales::trans_format("log10", scales::math_format(10^.x))
  #       ) +
  #       geom_hline(yintercept = EC_Lsmall, linetype = "dashed", size = .75) + annotate("text", x = res$Ltot %>% last(), y = as.numeric(EC_Lsmall), hjust = 0.75, vjust = 1.25, label = "EC['lim']", size = 5, parse = TRUE) + 
  #       labs(x = "Total ligand concentration (in nM)", y = "Concentration of tracer and complexes (in nM)") + theme(text = element_text(size=14)) +
  #       scale_color_manual(name = NULL, values = group.colors, labels = c("E", "EC", "C", "Ctot")) 
  #   } else { 
  #     pk <- res %>% select(Lfree, E, EC, Cfree, Ctot) %>% melt(id.vars = 'Lfree') %>% 
  #       ggplot(aes(x = Lfree, y = value, color = variable)) + geom_line(alpha=1, size=1.25) +
  #       annotation_logticks(sides = "b") + scale_x_continuous(trans="log10", #limits=c(1e-1,NA),
  #                                                             breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                                                             labels = scales::trans_format("log10", scales::math_format(10^.x))
  #       ) +
  #       geom_hline(yintercept = EC_Lsmall, linetype = "dashed", size = .75) + annotate("text", x = res$Lfree %>% last(), y = as.numeric(EC_Lsmall), hjust = 0.75, vjust = 1.25, label = "EC['lim']", size = 5, parse = TRUE) + 
  #       labs(x = "Free ligand concentration (estimated) (in nM)", y = "Concentration of tracer and complexes (in nM)") + theme(text = element_text(size=14)) +
  #       scale_color_manual(name = NULL, values = group.colors, labels = c("E", "EC", "C", "Ctot")) 
  #   }
  #   
  #   print(pk)
  # })
  
  # plotly
  output$plot_EC_E_E_Ctot <- renderPlotly({
    res = res()
    EC_Lsmall <- res$EC_Lsmall[1]
    
    if (res$Display_Ltot[1] == 1){
      pk <- res %>% select(Ltot, Efree, EC, Cfree, Ctot) %>% dplyr::rename(Cfree = Efree, Trfree = Cfree, CTr = EC, Trtot = Ctot) %>% melt(id.vars = 'Ltot') %>% 
        ggplot(aes(x = Ltot, y = value, color = variable)) + geom_line(alpha=1, size=1.25) +
        scale_x_continuous(trans="log10", breaks = scales::trans_breaks("log10", function(x) 10^x)) +
        geom_hline(yintercept = EC_Lsmall, linetype = "dashed", size = .75) + annotate("text", x = res$Ltot %>% last()*0.5, y = as.numeric(EC_Lsmall)*1.1, hjust = 0.75, vjust = 0.25, label = "EClim", size = 4, parse = TRUE) + 
        labs(x = "Total ligand concentration (in nM)", y = "Concentration of tracer and complexes (in nM)") +
        #scale_color_manual(name = NULL, values = group.colors, labels = c("E", "EC", "C", "Ctot")) 
        scale_color_manual(name = NULL, values = group.colors, labels = c("C", "C-Tr", "Tr", "Trtot")) 
    } else { 
      pk <- res %>% select(Lfree, Efree, EC, Cfree, Ctot) %>% dplyr::rename(Cfree = Efree, Trfree = Cfree, CTr = EC, Trtot = Ctot) %>% melt(id.vars = 'Lfree') %>% 
        ggplot(aes(x = Lfree, y = value, color = variable)) + geom_line(alpha=1, size=1.25) +
        scale_x_continuous(trans="log10", breaks = scales::trans_breaks("log10", function(x) 10^x)) +
        geom_hline(yintercept = EC_Lsmall, linetype = "dashed", size = .75) + annotate("text", x = res$Ltot %>% last()*0.5, y = as.numeric(EC_Lsmall)*1.1, hjust = 0.75, vjust = 0.25, label = "EClim", size = 4, parse = TRUE) + 
        labs(x = "Free ligand concentration (estimated) (in nM)", y = "Concentration of tracer and complexes (in nM)") + 
        #scale_color_manual(name = NULL, values = group.colors, labels = c("E", "EC", "C", "Ctot")) 
        scale_color_manual(name = NULL, values = group.colors, labels = c("C", "C-Tr", "Tr", "Trtot")) 
    }
    
    pk <- ggplotly(pk) # %>% layout(xaxis = list(exponentformat = "power"))
    print(pk)
  })
  
  output$param <- renderTable({
    if(input$Cooperativity_choice == 'alpha' && input$Native_interaction == 'native') {
      alpha_comp = input$Cooperativity 
      KP1_comp = ifelse(input$KT_choice == "KT1", input$KP1, input$KP1*alpha_comp)
      #KL = input$KP1*alpha_comp/input$KPE # input$KP1*input$KE1/(input$KE1/alpha_comp)/input$KPE
      KL = KP1_comp*alpha_comp/input$KPE # input$KP1*input$KE1/(input$KE1/alpha_comp)/input$KPE
      beta_comp = alpha_comp*input$KPE/input$KE1
      #gamma_comp = alpha_comp*input$KPE/input$KP1
      gamma_comp = alpha_comp*input$KPE/KP1_comp
    } else if(input$Cooperativity_choice == 'beta' && input$Native_interaction == 'native') {
      beta_comp = input$Cooperativity
      alpha_comp = beta_comp*input$KE1/input$KPE
      KP1_comp = ifelse(input$KT_choice == "KT1", input$KP1, input$KP1*alpha_comp)
      #KL = input$KP1/beta_comp
      KL = KP1_comp/beta_comp
      gamma_comp = input$KE1/KL
    } else if(input$Cooperativity_choice == 'gamma' && input$Native_interaction == 'native') {
      gamma_comp = input$Cooperativity 
      #alpha_comp = gamma_comp*input$KP1/input$KPE
      KP1_comp = NA
      alpha_comp = gamma_comp*KP1_comp/input$KPE
      #KL = input$KE1*input$KP1/alpha_comp/input$KPE #input$KE1*input$KP1/gamma_comp
      KL = input$KE1*KP1_comp/alpha_comp/input$KPE #input$KE1*input$KP1/gamma_comp
      #beta_comp = input$KP1/KL
      beta_comp = KP1_comp/KL
    } else {
      KL = 0
      alpha_comp = input$alpha # ifelse(input$alpha > 0, input$alpha, NULL) #
      KP1_comp = ifelse(input$KT_choice == "KT1", input$KP1, input$KP1*alpha_comp)
      beta_comp = NA
      gamma_comp = NA
    }
    
    data <- data.frame(Value = c(#input$KP1,
                                 KP1_comp,
                                 input$KE1, 
                                 ifelse(alpha_comp > 0 && is.finite(alpha_comp), alpha_comp, "Not defined"), 
                                 #ifelse(alpha_comp > 0 && is.finite(alpha_comp), input$KP1/alpha_comp, "Not defined"), 
                                 ifelse(alpha_comp > 0 && is.finite(alpha_comp), KP1_comp/alpha_comp, "Not defined"), 
                                 ifelse(alpha_comp > 0 && is.finite(alpha_comp), input$KE1/alpha_comp, "Not defined"), 
                                 ifelse(input$Native_interaction == 'none', 0, ifelse(input$KPE >= 0 && input$Cooperativity > 0, input$KPE, "Not defined")), 
                                 ifelse(input$Native_interaction == 'none', "Not defined", ifelse(KL >= 0 && input$Cooperativity >= 0 && input$KPE > 0, KL, "Not defined")), 
                                 ifelse(input$Native_interaction == 'none', 0, ifelse(beta_comp >= 0 && input$KPE >= 0, beta_comp, "Not defined")), 
                                 ifelse(input$Native_interaction == 'none', 0, ifelse(gamma_comp >= 0 && input$KPE >= 0, gamma_comp, "Not defined"))),
               row.names = c("<em>K</em><sub>T,1</sub>", "<em>K</em><sub>C,1</sub>", "α", "<em>K</em><sub>T,2</sub>", "<em>K</em><sub>C,2</sub>", "<em>K</em><sub>TC</sub>", "<em>K</em><sub>L</sub>", "β", "γ"))}, rownames = T, digits = 5, sanitize.text.function = function(x) x)
    
}
