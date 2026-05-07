

# Summarize predictions.
pred_summary <- function(mod_object,pars){
  pred_out <- list()
  for(i in 1:length(pars)){
    A <- extract(mod_object,pars=pars[i], permuted = TRUE)
    tmp_mean <- apply(A[[pars[i]]],c(2,3),mean) %>% as.data.frame() %>% 
      mutate(idx = 1:nrow(.))
    tmp_SD   <-  apply(A[[pars[i]]],c(2,3),sd) %>% as.data.frame() %>% 
      mutate(idx = 1:nrow(.))
    tmp_q025 <- apply(A[[pars[i]]],c(2,3),quantile,probs=0.025) %>% as.data.frame() %>% 
      mutate(idx = 1:nrow(.))
    tmp_q975 <- apply(A[[pars[i]]],c(2,3),quantile,probs=0.975) %>% as.data.frame() %>% 
      mutate(idx = 1:nrow(.))
    
    COL <- ncol(tmp_mean) - 1
    colnames(tmp_mean)[1:COL] <- paste0("sp_",1:COL)
    colnames(tmp_SD)   <- colnames(tmp_mean)
    colnames(tmp_q025) <- colnames(tmp_mean)
    colnames(tmp_q975) <- colnames(tmp_mean)
    
    tmp_mean2 <- pivot_longer(tmp_mean,cols = -idx, names_to= "species",
                              values_to=paste0("Mean")) %>% mutate(param = pars[i])
    tmp_SD2 <- pivot_longer(tmp_SD,cols = -idx, names_to= "species",
                            values_to=paste0("SD"))
    tmp_q025a <- pivot_longer(tmp_q025,cols = -idx, names_to= "species",
                              values_to=paste0("q025"))
    tmp_q975a <- pivot_longer(tmp_q975,cols = -idx, names_to= "species",
                              values_to=paste0("q975"))
    tmp_long <- left_join(tmp_mean2,tmp_SD2) %>% 
      left_join(.,tmp_q025a) %>% 
      left_join(.,tmp_q975a)
    
    pred_out[[pars[i]]] <- 
      list(Mean = tmp_mean,
           SD  = tmp_SD,
           q025 = tmp_q025,
           q975 = tmp_q975,
           Long = tmp_long
      )
  }
  return(pred_out)
}