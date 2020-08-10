function [pCOD1, pCOD2, pCOD12, ss_res1, ss_res2, ss_res12, ss_total] = ABCD_pCOD_2groups ...
    (y_pred1, y_pred2, y_test1, y_test2, y_train1, y_train2)

    ss_res1 = sum((y_pred1 - y_test1).^2, 1) ./ length(y_test1);
    ss_res2 = sum((y_pred2 - y_test2).^2, 1) ./ length(y_test2);
    ss_res12 = sum( ([y_pred1; y_pred2] - [y_test1; y_test2]).^2, 1 ) ...
        ./ length([y_test1; y_test2]);
    ss_total = sum( ([y_train1; y_train2] - mean([y_train1; y_train2],1)).^2, 1 ) ...
        ./ length([y_train1; y_train2]);
    
    pCOD1 = bsxfun(@minus, 1, ss_res1./ss_total);
    pCOD2 = bsxfun(@minus, 1, ss_res2./ss_total);
    pCOD12 = bsxfun(@minus, 1, ss_res12./ss_total);
end