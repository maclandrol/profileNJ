function r=compar(x, y)
r=zeros(size(x));
for i=1:size(x,1)
    r(i, :)= x(i, :)==y(i);
    r(i, r(i,:)~=1)=x(i, r(i,:)~=1)/y(i);
    
end

end