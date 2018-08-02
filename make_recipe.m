#Denominator reduction
demred:=proc(x, sequence := [], useorder := false) local j1,j2,thedegree_term:
print("reduced"):
j1:=indets(x):
if j1=[] then return x else j1:=j1[1] end if:

if useorder then 
    if sequence = [] then return x end if:
    j1 := sequence[1]:
end if:

thedegree_term := degree(x,j1):

#may happen if we choose a bad variable
if thedegree_term= 0 then 
    return demred(x,sequence[2..],true):
end if:

if thedegree_term>2 then error "%1 is not quadratic in %2",x,j1 end if:
if thedegree_term =1 then 
    if useorder then 
        return demred(factor(coeff(x,j1)),sequence[2..],true):
    else 
        return demred(factor(coeff(x,j1))):
    end if:
end if:
print(nops(x)):
j2:=map(proc(y) if depends(y,j1) then y else NULL end if end proc,[op(x)])[1]:
if not (divide(x,j2)) then return x end if:
if degree(j2,j1) > 1 then return x end if: 
    #print("resultant",nops(j2),nops(x/j2),nops(x),j1):
j2:=factor(resultant(j2,x/j2,j1)):
j1:=indets(j2):

if j1=[] then return j2
else 
    if useorder then 
        if sequence[2..] <> [] then
            j1 := sequence[2]
        else
            j1 := j1[1]
        end if:
    else
        j1 := j1[1]
    end if:
end if:

if max(map(y->degree(y,j1),[op(j2)]))>1 then x
else 
    if useorder then
        demred(j2,sequence[2..],true) :
    else 
        return demred(j2):
    end if:
end if 
end proc:


with (LinearAlgebra):
read(fiveinv):
nthinv := fiveinv:
print(nops(indets(nthinv))):
read(elim_seq):

print(nops(nthinv)):
nthinv := demred(factor(nthinv),elim_seq,true):

#homogenous
#nthinv := subs(indets(nthinv)[-1]=1,nthinv):


print(nops(indets(nthinv))):
primes := {3}:
print(nops(expand(nthinv))):
ti := time():

##Enter a file location inside the quotation marks.
recipename := cat("",convert(graphname,string)):
print(recipename):
fd := fopen(recipename, WRITE):
################################
  btopoly := table():
  btopol_reverse := table():
  modded := nthinv:
  factinv := factors(modded)[2]:
  #print(factors(modded)[1]):
  #print(factinv):
  count := 1:

  list_sixinv := []:
  sixinv := 1:
  for fact in factinv do
    polyfact:= fact[1]:
    exponent := fact[2]:
    polyfact := expand(polyfact^exponent):
    btopoly[polyfact] := cat(b_,count):
    btopol_reverse[cat(b_,count)] := polyfact:
    sixinv := sixinv * cat(b_,count):
    list_sixinv := [op(list_sixinv), cat(b_,count)]:
    count := count + 1:
   end do:

   sixinv := sixinv * factors(modded)[1]:

sixinv := sixinv^(5-1):
print(sixinv):
all_indets := indets(expand(nthinv)):
most_terms := 0:
loopcount := 1:

for j in elim_seq do
    var_to_extract := j:
    if(not has(all_indets,var_to_extract)) then next end if:
    all_indets := subs(var_to_extract = NULL,all_indets):
    print(op(list_sixinv), "list sixinv size"): 
    temp_list_sixinv := []:
    toprint_list := []:
    for var in list_sixinv do
        new_poly1 := coeff(btopol_reverse[var],var_to_extract):
        new_poly2 := subs(var_to_extract=0,btopol_reverse[var]):

        tosub1 := 1:
        tosub2 := 1:

        if new_poly1 = 0 then
            tosub1 := 0:
        end if:

        if new_poly2 = 0 then
            tosub2 := 0:
        end if:

       for fact in factors(new_poly1)[2] do
           if assigned(btopoly[fact[1]]) then 
               tosub1 := tosub1 * btopoly[fact[1]]^fact[2]:
           else 
               btopoly[fact[1]] := cat(b_,count):
               btopol_reverse[cat(b_,count)] := fact[1]:
               tosub1 := tosub1 * btopoly[fact[1]]^fact[2]:
               count := count + 1:
           end if:

           if not member(btopoly[fact[1]], temp_list_sixinv) then
               temp_list_sixinv := [op(temp_list_sixinv), btopoly[fact[1]]]:
           end if:
       end do:

       tosub1 := tosub1 * factors(new_poly1)[1]:

       for fact in factors(new_poly2)[2] do
           if assigned(btopoly[fact[1]]) then 
               tosub2 := tosub2 * btopoly[fact[1]]^fact[2]:
           else 
               btopoly[fact[1]] := cat(b_,count):
               btopol_reverse[cat(b_,count)] := fact[1]:
               tosub2 := tosub2 * btopoly[fact[1]]^fact[2]:
               count := count + 1:
           end if:
            if not member(btopoly[fact[1]], temp_list_sixinv) then
               temp_list_sixinv := [op(temp_list_sixinv), btopoly[fact[1]]]:
           end if:
       end do:

       tosub2 := tosub2 * factors(new_poly2)[1]:
       print(var,tosub1,tosub2):
       toprint_list := [op(toprint_list), cat(convert(var,string),",",convert(tosub1,string),",",convert(tosub2,string))];

   end do:

   fprintf(fd,">%d\n",nops(list_sixinv)):

   for v in op(list_sixinv) do
       print(v):
       fprintf(fd, "-%s\n",convert(v,string)):
   end do:

   for v in op(temp_list_sixinv) do
       fprintf(fd, "+%s\n",convert(v,string)):
   end do:

   for str in toprint_list do
       fprintf(fd,"%s\n",str):
   end do:
   list_sixinv := temp_list_sixinv:
end do:
########################
fclose(fd):
quit:
