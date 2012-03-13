function makenewtonplot%(d, n, k, incr)
%     figure
%     for i = 1:k
%         m = i * incr;
%         subplot(k, 2, 2 * i - 1)
        fastnewton(30);
%         title(sprintf('$d=%i$', m), 'interpreter', 'latex')
%     end
    
%     subplot(k, 2, 2:2:(k * 2))
    
    % run the tester n times solving from 2x2 up to dxd matrices
%     plot(test_tdma(n, d))
%     xlabel('$d$', 'interpreter', 'latex')
%     ylabel('Time', 'interpreter', 'latex')
%     title(['$d_\mathrm{max}=' num2str(d) '$'], 'interpreter', 'latex')
%     axis tight
%     pdfsave('newtons_method.pdf')
end
