function flag = isPositiveInt(x)
if all(x>0 & mod(x,1)==0)
    flag = true;
else
    flag = false;
end
end