function [err1, err2, err12] = ABCD_errors_2groups(y_pred1, y_pred2, y_test1, y_test2, metric, y_train1, y_train2)

% [err1, err2, err12] = ABCD_errors_2groups(y_pred1, y_pred2, y_test1, y_test2, metric, y_train1, y_train2)
%
% Given the predicted scores `y_pred?`, the raw test scores `y_test?`, and the raw training scores 
% `y_train?`, calculate the error-based accuracy metrics. 
% `metric`: choose from: 'MSE', 'MSE_norm', 'MAE', 'MAE_norm'


err1 = CBIG_compute_prediction_acc_and_loss(y_pred1, y_test1, metric, y_train1);
err2 = CBIG_compute_prediction_acc_and_loss(y_pred2, y_test2, metric, y_train2);
err12 = CBIG_compute_prediction_acc_and_loss([y_pred1; y_pred2], [y_test1; y_test2], metric, [y_train1; y_train2]);
    
end