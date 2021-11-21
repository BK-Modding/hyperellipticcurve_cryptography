(* ::Package:: *)

(* This package was created as a part of diploma thesis: *)
(* Hyperelliptic curves and their application in cryptography *)

BeginPackage["HyperellipticCryptography`HyperellipticCryptography`"]
OrdDiv::usage="..."
MultiDiv::usage="..."
KeyGen::usage="..."
Encrypt::usage="..."
Decrypt::usage="..."

(* >> HyperellipticCurves.m *)

SetDirectory[NotebookDirectory[]]
Import["HyperellipticCurves.m"]


Begin["`Private`"]

(*set the private context*)

protected=Unprotect[OrdDiv,MultiDiv,KeyGen,Encrypt,Decrypt]



OrdDiv[ap_,bp_,K_,k_,h_,f_,x_,g_]:=(
aq=ap;
bq=bp;
test=0;
n=1;
For[i=1,test==0,i++,{{aq,bq}=Cantor[aq,bq,ap,bp,h,f,x,g],n++,If[{aq,bq}=={k[1],0},test=1]}];
n
)

MultiDiv[a_,b_,multi_,K_,k_,h_,f_,x_,g_]:=(
bin=IntegerDigits[multi,2];
cc={};
For[i=1,i<Length[bin]+1,i++,AppendTo[cc,0]];
list={};
r={};
For[i=1,i<Length[bin]+1,i++,If[bin[[i]]==1,{cc[[i]]=1,AppendTo[list,cc],AppendTo[r,Length[bin]-i],cc[[i]]=0}]];
l={};
For[i=1,i<Length[list]+1,i++,AppendTo[l,FromDigits[list[[i]],2]]];
dDiv={};
For[i=1,i<Length[r]+1,i++,{{aa,bb}={a,b},For[j=1,j<r[[i]]+1,j++,{{aa,bb}=DoubleDiv[aa,bb,h,f,x],{aa,bb}=RedDiv[aa,bb,2,h,f,x]}],AppendTo[dDiv,{aa,bb}]}];
{ad,bd}={k[1],0};
For[i=1,i<Length[dDiv]+1,i++,{ad,bd}=Cantor[ad,bd,dDiv[[i]][[1]],dDiv[[i]][[2]],h,f,x,2]];
{ad,bd}
)

KeyGen[ap_,bp_,ord_,K_,k_,h_,f_,x_,g_]:=(
key=RandomInteger[{1,ord}];
{aq,bq}=MultiDiv[ap,bp,key,K,k,h,f,x,g];
{key,{aq,bq}}
)

Encrypt[am_,bm_,ap_,bp_,ord_,aq_,bq_,K_,k_,h_,f_,x_,g_]:=(
n=OrdDiv[ap,bp,K,k,h,f,x,g];
a=RandomInteger[{1,ord}];
{ac1,bc1}=MultiDiv[ap,bp,a,K,k,h,f,x,g];
{aac2,bbc2}=MultiDiv[aq,bq,a,K,k,h,f,x,g];
{ac2,bc2}=Cantor[am,bm,aac2,bbc2,h,f,x,g];
{{ac1,bc1},{ac2,bc2}}
)

Decrypt[ac1_,bc1_,ac2_,bc2_,key_,K_,k_,h_,f_,x_,g_]:=(
{ac1op,bc1op}={ac1,PolynomialRemainder[-h-bc1,ac1,x]};
{akc1,bkc1}=MultiDiv[ac1op,bc1op,key,K,k,h,f,x,g];
{amFin,bmFin}=Cantor[ac2,bc2,akc1,bkc1,h,f,x,g];
{amFin,bmFin}
)


Protect[Evaluate[protected]]        (*restore protection of system symbols*)
End[]                               (*set the private context*)
EndPackage[]                        (*end the package context*)
