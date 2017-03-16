Ac=c(1/3,1/3,1/3,
     1/3,1/3,1/3,
     1/3,1/3,1/3);
Bc=c(1/4,1/4,1/4,1/4,0,0,0,0,
     1/6,1/6,1/6,1/6,1/6,1/6,0,0,
     1/8,1/8,1/8,1/8,1/8,1/8,1/8,1/8);
A=matrix(Ac,nrow=3,byrow=T);
B=matrix(Bc,nrow=3,byrow=T);
colnames(A)=c("D4","D6","D8");
rownames(A)=c("D4","D6","D8");
colnames(B)=c("1","2","3","4","5","6","7","8");
rownames(B)=c("D4","D6","D8");
PI=c(1/3,1/3,1/3);
O=c("1","6","3","5","2","7","3","5","2","4");





forward_procedure=function(O,A,B,PI)
{
  alpha=matrix(0,nrow=length(O),ncol=dim(B)[1]);
  for(t in 1:length(O))
  {
    for(j in 1:dim(B)[1])
    {
      if(t==1)
      {
        alpha[t,j]=PI[j]*B[j,O[t]];
        cat(paste("alpha[",t,",",j,"]","=",alpha[t,j],"\n",sep=""));
      }else{
        alpha[t,j]=alpha[t-1,]%*%A[,j]*B[j,O[t]];
        cat(paste("alpha[",t,",",j,"]","=",alpha[t,j],"\n",sep=""));
      }
    }
  }
  p=sum(alpha[length(O),]);
  return(p);
}

p=forward_procedure(O,A,B,PI);
cat(paste("p(O|mu)=",p,"\n",sep=""));




viterbi=function(O,A,B,PI)
{
  delta=matrix(0,nrow=length(O),ncol=dim(A)[1]);
  phi=matrix(0,nrow=length(O),ncol=dim(A)[1]);
  for(t in 1:length(O))
  {
    for(i in 1:dim(A)[1])
    {
      if(t==1)
      {
        delta[t,i]=PI[i]*B[i,O[t]];
        phi[t,i]=0;
        #        cat(paste("delta[",t,",",i,"]=",delta[t,i],"\n",sep=""));
        #        cat(paste("phi[",t,",",i,"]=",phi[t,i],"\n",sep=""));
      }else{
        delta_a=delta[t-1,]*A[,i];
        delta[t,i]=max(delta_a)*B[i,O[t]];
        phi[t,i]=which.max(delta_a);
        #        cat(paste("delta[",t,",",i,"]=",delta[t,i],"\n",sep=""));
        #        cat(paste("phi[",t,",",i,"]=",phi[t,i],"\n",sep=""));
      }     
    }
  }
  colnames(delta)=colnames(A);
  rownames(delta)=1:length(O);
  colnames(phi)=colnames(A);
  rownames(phi)=1:length(O);
  cat("\nThe calculation is complete.\n")
  cat("\ndelta_matrix:\n");
  print(delta);
  cat("\nphi_matrix:\n");
  print(phi);
  q=rep(0,length(O));
  cat("\nStart backtracking:\n");
  for(t in seq(from=length(O),to=1,by=-1))
  {
    if(t==length(O))
    {
      q[t]=which.max(delta[length(O),]);
      cat(paste(colnames(A)[q[t]],"<-",sep=""));
    }else{
      q[t]=phi[t+1,q[t+1]];
      if(t!=1)
        cat(paste(colnames(A)[q[t]],"<-",sep=""));
      if(t==1)
        cat(paste(colnames(A)[q[t]],"\n",sep=""));
    }
  }
  return(colnames(A)[q]);
}

Q=viterbi(O,A,B,PI);
cat(c(Q,"\n"));








O=unlist(strsplit("ÑÐ¾¿ÉúÎïºÜÓÐÒâË¼Ëû´óÑ§Ê±´úÊÇÑÐ¾¿ÉúÎïµÄÉúÎï×¨ÒµÊÇËûµÄÊ×Ñ¡Ä¿±ê",split=""));
Q=c("B","E","B","E","S","B","M","E",
    "S","B","E","B","E","S","B","E","B","E","S",
    "B","E","B","E","S","S","S","B","M","M","E");




Delta=function(Q,S)
{
  delta=matrix(0,nrow=length(Q),ncol=length(S));
  for(i in 1:length(Q))
  {
    for(j in 1:length(S))
    {
      delta[i,j]=ifelse(Q[i]==S[j],1,0);
    }
  }
  return(delta);
}

MLE=function(O,Q)
{
  S=unique(Q);
  K=unique(O);
  A=matrix(0,nrow=length(S),ncol=length(S));
  B=matrix(0,nrow=length(S),ncol=length(K));
  PI=rep(0,length(S));
  T=length(O);
  delta_q_s=Delta(Q,S);
  delta_O_v=Delta(O,K);
  for(i in 1:length(S))
  {
    PI[i]=delta_q_s[1,i];
    for(j in 1:length(S))
    {
      delta_qt_si_a=delta_q_s[1:(T-1),i];
      delta_qt1_sj_a=delta_q_s[2:T,j];
      delta_si_sj_a=delta_qt_si_a*delta_qt1_sj_a;
      A[i,j]=sum(delta_si_sj_a)/sum(delta_qt_si_a);
      
    }
    for(k in 1:length(K))
    {
      delta_qt_si_b=delta_q_s[,i];
      delta_Ot_vk_b=delta_O_v[,k];
      delta_si_vk_b=delta_qt_si_b*delta_Ot_vk_b;
      B[i,k]=sum(delta_si_vk_b)/sum(delta_qt_si_b);
    }
  }
  names(PI)=S;
  colnames(A)=S;
  rownames(A)=S;
  colnames(B)=K;
  rownames(B)=S;
  result=list();
  result$PI=PI;
  result$A=A;
  result$B=B;
  return(result);
}

mu=MLE(O,Q);
print(mu);


