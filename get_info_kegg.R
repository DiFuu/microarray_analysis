
library(KEGGREST)

############################################################
### get_info kegg : get the description of a series of pathway ids passed as argument
## pathways_ids: array of pathway identifiers
## the function returns, in a single string and for all pathways, the pathway identifier
##  followed by the information
############################################################

myKeggGet<-function(keggId)
{
  out <- tryCatch(
    {
      #KEGGREST function that returns for a pathway ID, its description (NAME attribute):
      keggGet(keggId)[[1]]$NAME
    },
    error=function(cond) {
      message(paste("Identifier",keggId, "not found"))
      return('')
    },
    finally={}
  )  
  return(out)

}



get_info_kegg<-function(pathway_ids)
{
  
  info_kegg<-c()
  
  if (any(is.na(pathway_ids))==FALSE) #If there are pathway identifiers
  {
    for(current_pathway in pathway_ids) #We go through all the pathway ids
    {
      print(paste('Getting pathway information for',current_pathway))
     
      kegg_name <- myKeggGet(paste('hsa',current_pathway,sep=''))
      
      #We concatenate the result in the form <pathway_id>:info_pathway
      info_kegg<-cbind(info_kegg,paste(current_pathway,kegg_name,sep=':'))
    }
  }
  
  return(paste(info_kegg,collapse=' // '))
  
}
