function mysaveas(fig,fn)
if isMatlab,
  saveas(fig,sprintf('%s.pdf',fn));
  saveas(fig,sprintf('%s.jpg',fn));
end;
