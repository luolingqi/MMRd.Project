
buildLmeModel <- function(x,y,z,data) {# x: response; y: fixed_effects; z: random_effects - cancerType
  # make categorical
  data[,y[2]] <- to_factor(data[,y[2]]) # MMRd_Type
  print(data[,y[2]])
  data[,y[1]] <- to_factor(data[,y[1]]) # Cancer_Type_NEW
  print(data[,y[1]])
  ###############################################################################
  # simple model without any interaction terms
  ###############################################################################
  formula <- reformulate(termlabels = c(y, paste0("(1+",y[1],"|",z,")")), 
                         response = x)
  #print(paste0(format(formula), collapse = ""))
  
  fit0 <- lmer(formula,
               data = data, REML = FALSE)
  #print(get(paste0(response,".model")))
  sink(file = file.path(getwd(),"results_model_cmp",paste0("summary.",x,".lmer.model.W.fixed.",paste0(y,collapse = "."),
                                                           ".random_", paste0(z,collapse = "."),"_Sample_Type.txt")),
       append = F)
  print(fit0)
  print(coef(fit0))
  sink()
  
  # plot the model (fixed effect - MMRd type)
  p <- plot_model(fit0, sort.est = TRUE,show.values = TRUE, value.offset = .3,
                  title = paste0("linear mixed effect model for ",x,"\n",gsub("\\(","\n\\(", paste0(format(formula), collapse = "")))) +
    theme(plot.title = element_text(size = 10)) + theme(axis.text.y = element_text(size = 10))
  print(p)
  # # plot the model (interaction)
  # p <- plot_model(fit0, type = "pred", terms = c(y[1],y[2])) + theme(axis.text.x = element_text(angle = 45,hjust = 1))
  # print(p)
  
  
  ###############################################################################  
  # model with interaction term between cancer type and MMRd type
  ###############################################################################
  #formula <- reformulate(termlabels = c(y, paste0(y[2],"*",y[1]),  paste0("(1+",y[1],"|",z,")")), response = x)
  formula <- reformulate(termlabels = c(paste0(y[2],"*",y[1]), y[3:length(y)], paste0("(1+",y[1],"|",z,")")), response = x)
  #formula <- reformulate(termlabels = c(paste0(y[2],"*",y[1]), y[3:length(y)]), response = x)
  #print(paste0(format(formula), collapse = ""))
  
  fit1 <- lmer(formula,
               data = data, REML = FALSE)
  #print(get(paste0(response,".model")))
  sink(file = file.path(getwd(),"results_model_cmp",paste0("summary.",x,".lmer.model.W.fixed.",paste0(y,collapse = "."),".w.interaction",
                                                           ".random_", paste0(z,collapse = "."),"_Sample_Type.txt")),
       append = F)
  print(fit1)
  print(coef(fit1))
  sink()
  
  # # plot the model (fixed effect - MMRd type)
  # p <- plot_model(fit1, sort.est = TRUE,show.values = TRUE, value.offset = .3,
  #                 title = paste0("linear mixed effect model for ",x,"\n",gsub("\\(","\n\\(", paste0(format(formula), collapse = "")))) +
  #   theme(plot.title = element_text(size = 8)) + theme(axis.text.x = element_text(size = 5))
  # print(p)
  # plot the model (interaction)
  #data[,y[2]]
  p <- plot_model(fit1, type = "eff", terms = c(y[1],y[2]),
  #p <- plot_model(fit1, type = "int",
                  title = paste0("linear mixed effect model for ",x,"\n",gsub("\\(","\n\\(", paste0(format(formula), collapse = "")))) + 
    theme(axis.text.x = element_text(angle = 45,hjust = 1)) + theme(plot.title = element_text(size = 8))
  print(p)
  
  
  ############################################################################### 
  # model w/o random effects
  ###############################################################################
  #formula <- reformulate(termlabels = c(y,paste0(y[2],"*",y[1])),response = x)
  formula <- reformulate(termlabels = c(paste0(y[2],"*",y[1]), y[3:length(y)]),response = x)
  #print(paste0(format(formula), collapse = ""))
  
  fit2 <- lm(formula, 
             data = data)
  
  #print(get(paste0(response,".model")))
  sink(file = file.path(getwd(),"results_model_cmp",paste0("summary.",x,".lmer.model.W.fixed.",paste0(y,collapse = "."),
                                                           ".WO.random_", paste0(z,collapse = "."),"_Sample_Type.txt")),
       append = F)
  print(fit2)
  print(coef(fit2))
  sink()
 
  #plot the model (fixed effect - MMRd)
  # p <- plot_model(fit2,sort.est = TRUE,show.values = TRUE, value.offset = .3,
  #                 title = paste0("linear mixed effect model for ",x,"\n",gsub("\\(","\n\\(", paste0(format(formula), collapse = "")))))) +
  #   theme(plot.title = element_text(size = 8)) + theme(axis.text.x = element_text(size = 5))
  # print(p)
  # #plot the model (interaction)
  # p <- plot_model(fit2, type = "pred", terms = c(y[1],y[2])) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # print(p)
  
  
  
  # print anova test for the comparison of the above 2 models W/WO fixed effect
  sink(file = file.path(getwd(),"results_model_cmp",paste0("summary.",x,".lmer.model.W.fixed.",paste0(y,collapse = "."),
                                                           ".W.vs.WO.random_",paste0(z,collapse = "."),"_Sample_Type_anova.txt")),
       append = F)
  print(anova(fit1,fit2))
  sink()

  
  return(list(fit1,fit2))
}