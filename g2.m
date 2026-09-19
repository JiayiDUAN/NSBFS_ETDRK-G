function g2=g2(phi)
g2=(gradx(phi).^2+gradx_m(phi).^2+gradyNeu(phi).^2+grady_mNeu(phi).^2)/2;
end


