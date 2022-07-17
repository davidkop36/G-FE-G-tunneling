function I= ExtractMinimumInDirection (x,dir)

[~,I] = min(x,[],dir);
I = squeeze(I);
end