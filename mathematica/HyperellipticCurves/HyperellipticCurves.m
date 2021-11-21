(* ::Package:: *)

(* This package was created as a part of diploma thesis: *)
(* Hyperelliptic curves and their application in cryptography *)

BeginPackage["HyperellipticCurves`HyperellipticCurves`",{"FiniteFields`FiniteFields`","Algebra`PolynomialExtendedGCD`"}]

SetField::usage="..."
FindPoints::usage="..."
FindPointsIm::usage="..."
NumberOfPoints::usage="..."
OrdOrdinary::usage="..."
OrdSpecial::usage="..."
OrdInfinity::usage="..."
NormG::usage="..."
AddDiv::usage="..."
DoubleDiv::usage="..."
RedDiv::usage="..."
JacobianOrd::usage="..."
DivToPol::usage="..."
PolToDiv::usage="..."
Cantor::usage="..."





<<FiniteFields`FiniteFields`


Begin["`Private`"]

(*set the private context*)

protected=Unprotect[SetField,FindPoints,FindPointsIm,NumberOfPoints,OrdOrdinary,OrdSpecial,OrdInfinity,NormG,AddDiv,DoubleDiv,RedDiv,JacobianOrd,DivToPol,PolToDiv,Cantor]




SetField[p_,n_,K_,k_,rep_]:=(
K=GF[p,n];
r=Total[ToCharacterCode[rep]];
If[r==318,SetFieldFormat[K,FormatType->FunctionOfCode[k]],
If[r==331,SetFieldFormat[K,FormatType->FunctionOfCoefficients[k]],
If[r==313,SetFieldFormat[K,FormatType->Subscripted],Print["Unexpected input."]]]]
)

FindPoints[K_,k_,h_,f_,x_]:=(
Points={};
For[i=0,i<FieldSize[K],i++,For[j=0,j<FieldSize[K],j++,
fa=f/.x->k[i];
ha=h/.x->k[i];
If[fa==k[j]^2+k[j]*ha,Points=Union[Points,{{i,j}}]]]];
Points
)

FindPointsIm[K_,k_,h_,f_,x_]:=(
points={};
For[i=0,i<FieldSize[K],i++,For[j=0,j<FieldSize[K],j++,
fa=f/.x->k[i];
ha=h/.x->k[i];
If[fa==k[j]^2+k[j]*ha,{AppendTo[points,{k[i],k[j]}], AppendTo[points,{k[i],-k[j]-ha}],Break[]}]]];
Points=Union[points] ;
Points
)

NumberOfPoints[K_,k_,h_,f_,x_]:=(
1+Length[FindPoints[K,k,h,f,x]]
)

OrdOrdinary[h_,f_,P_,a_,b_,x_]:=(
pra=0;ra=0;prb=0;rb=0;s=0;pr=0;
For[i=1,pra==0,i++,{pra=PolynomialRemainder[a,(x-P[[1]])^i,x],ra++}];
For[i=1,prb==0,i++,{prb=PolynomialRemainder[b,(x-P[[1]])^i,x],rb++}];
If[ra>rb,r=rb-1,r=ra-1];
aa=PolynomialQuotient[a,(x-P[[1]])^r,x];
bb=PolynomialQuotient[b,(x-P[[1]])^r,x];
NormG0=aa^2+aa*bb*h-bb^2*f;
For[i=1,pr==0,i++,{pr=PolynomialRemainder[NormG0,(x-P[[1]])^i,x],s++}];
OrdG=r+s-1;
OrdG
)

OrdSpecial[h_,f_,P_,a_,b_,x_]:=(
pra=0;ra=0;prb=0;rb=0;s=0;pr=0;
For[i=1,pra==0,i++,{pra=PolynomialRemainder[a,(x-P[[1]])^i,x],ra++}];
For[i=1,prb==0,i++,{prb=PolynomialRemainder[b,(x-P[[1]])^i,x],rb++}];
If[ra>rb,r=rb-1,r=ra-1];
aa=PolynomialQuotient[a,(x-P[[1]])^r,x];
bb=PolynomialQuotient[b,(x-P[[1]])^r,x];
NormG0=aa^2+aa*bb*h-bb^2*f;
For[i=1,pr==0,i++,{pr=PolynomialRemainder[NormG0,(x-P[[1]])^i,x],s++}];
OrdG=2*r+s-1;
OrdG
)

OrdInfinity[g_,a_,b_,x_]:=(
OrdG=-Max[2*Exponent[a,x],2*Exponent[b,x]+2*g+1];
OrdG
)

NormG[a_,b_,h_,f_,x_]:=(
aa=a;
bb=b;
NG=Expand[aa^2+aa*bb*h-bb^2*f];
NG
)

AddDiv[a1_,b1_,a2_,b2_,h_,f_,x_]:=(
{d1,{e1,e2}}=PolynomialExtendedGCD[a1,a2,x];
{d,{c1,c2}}=PolynomialExtendedGCD[d1,(b1+b2+h),x];
s1=Expand[c1*e1];
s2=Expand[c1*e2];
s3=c2;
a=PolynomialQuotient[a1*a2,d^2,x];
bb=PolynomialQuotient[s1*a1*b2+s2*a2*b1+s3*(b1*b2+f),d,x];
b=PolynomialRemainder[bb,a,x];
{a,b}
)

DoubleDiv[a1_,b1_,h_,f_,x_]:=(
{d,{s1,s3}}=PolynomialExtendedGCD[a1,(b1+b1+h),x];
a=PolynomialQuotient[a1*a1,d^2,x];
b=PolynomialRemainder[PolynomialQuotient[s1*a1*b1+s3*(b1*b1+f),d,x],a,x];
{a,b}
)

RedDiv[a_,b_,g_,h_,f_,x_]:=(
aa=a;
bb=b;
For[i=0,Exponent[aa,x]>g,i++,{a1=PolynomialQuotient[f-bb*h-bb^2,aa,x],b1=PolynomialRemainder[-h-bb,a1,x],aa=a1,bb=b1}];
ex=Exponent[aa,x];
inv=1/Coefficient[aa,x,ex];
aa=Expand[inv*aa];
{aa,bb}
)

JacobianOrd[K_,M1_,M2_,n_]:=(
p=Characteristic[K];
a1=M1-1-p;
a2=(M2-1-p^2+a1^2)/2;
Res=Solve[x^2+a1*x+(a2-2*p)==0,x];
gamma1=Res[[1]][[1]][[2]];
gamma2=Res[[2]][[1]][[2]];
Res1=Solve[x^2-gamma1*x+p==0,x];
Res2=Solve[x^2-gamma2*x+p==0,x];
alpha1=Res1[[1]][[1]][[2]];
alpha2=Res2[[1]][[1]][[2]];
Nn=Abs[1-alpha1^n]^2*Abs[1-alpha2^n]^2;
JacobOrd=Expand[Nn/.n->ExtensionDegree[K]];
JacobOrd
)

DivToPol[K_,k_,h_,f_,x_,Div_]:=(
l=Length[Div];
ad=Div[[All,2]];
aa=k[1];
For[i=1,i<Length[ad]+1,i++,aa=aa*(x-ad[[i]])^(Div[[i]][[1]])];
Pola=Expand[aa];
ap=Apart[1/(aa),x];
app={};
For[i=1,i<Length[ap]+1,i++,app=AppendTo[app,ap[[i]]]];
pow=Div[[All,1]];
m=Max[pow];
For[q=1,q<m+1,q++,For[j=1,j<Length[app]+1,j++,For[i=1,i<Length[app]+1,i++,If[i!=j,If[PolynomialRemainder[Denominator[app[[j]]],Denominator[app[[i]]],x]==0,{tog=Together[app[[i]]+app[[j]]],app=ReplacePart[app,j->tog],app=Delete[app,i],i=i-1}]]]]];
nDiv={};
For[i=1,i<Length[Div]+1,i++,{{p,q}={x-Div[[i]][[2]],Div[[i]][[3]]},{pn,qn}={p,q},If[Div[[i]][[1]]>1,For[j=1,j<Div[[i]][[1]],j++,{pn,qn}=AddDiv[p,q,pn,qn,h,f,x]]],AppendTo[nDiv,{pn,qn}]}];
den={};
For[i=1,i<Length[app]+1,i++,AppendTo[den,Denominator[app[[i]]]]];
den=Expand[den];
a2={};
For[i=1,i<Length[den]+1,i++,For[j=1,j<Length[nDiv]+1,j++,If[den[[i]]==nDiv[[j]][[1]],AppendTo[a2,nDiv[[j]][[2]]*app[[i]]]]]];
ssuma={};
For[i=1,i<Length[a2]+1,i++,AppendTo[ssuma,a2[[i]]*aa]];
ssuma=Expand[ssuma];
suma=0;
For[i=1,i<Length[ssuma]+1,i++,suma=Expand[suma+ssuma[[i]]]];
b=PolynomialRemainder[suma,Pola,x];
{Pola,b}
)

PolToDiv[K_,k_,x_,a_,b_]:=(
elem=Table[k[i],{i,0,FieldSize[K]-1}];
aa=a;
xcor={};
For[i=1,Length[xcor]<Exponent[a,x],i++,{ae=aa/.x->elem,For[j=1,j<Length[ae]+1,j++,If[ae[[j]]==0,{AppendTo[xcor,k[j-1]],aa=PolynomialQuotient[aa,(x-k[j-1]),x]}]],If[xcor=={},{Print["Points of this divisor are not in the basic field"],Break[]}]}];
Div={};
For[i=1,i<Length[xcor]+1,i++,AppendTo[Div,{xcor[[i]],b/.x->xcor[[i]]}]];
For[i=1,i<Length[Div]+1,i++,{n=1;For[j=1,j<Length[Div]+1,j++,If[i!=j,If[Div[[i]][[1]]==Div[[j]][[1]],{n=n+1,Div=Delete[Div,j],j=j-1}]]],Div=ReplacePart[Div,i->{n,Div[[i]][[1]],Div[[i]][[2]]}]}];
For[i=1,i<Length[Div]+1,i++,If[Length[Div[[i]]]<3,Div=ReplacePart[Div,i->{1,Div[[i]][[1]],Div[[i]][[2]]}]]];
Div
)

Cantor[a1_,b1_,a2_,b2_,h_,f_,x_,g_]:=(
{a,b}=AddDiv[a1,b1,a2,b2,h,f,x];
{aa,bb}=RedDiv[a,b,g,h,f,x];
{aa,bb}
)


Protect[Evaluate[protected]]        (*restore protection of system symbols*)
End[]                               (*set the private context*)
EndPackage[]                        (*end the package context*)





