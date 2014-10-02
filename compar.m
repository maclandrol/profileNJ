function r=compar(x, y)
r=zeros(size(x));
for i=1:size(x,1)
    if(x(i)==y(i))
        r(i)=1;
    else
        r(i)=x(i)/y(i);
    end
        
end

end