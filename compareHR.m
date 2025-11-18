function [HR, index] = compareHR(refHR,rates, times)
    % Return the HR closest to refHR and the index of it from vector.
    % In case of tie, go base off time to run the code.

    % find minimum difference, smallest difference will be the HR closest
    % to the refHR
    diffs = abs(rates-refHR);
    [minDiff, index1] = min(diffs); %the index will be used to determine method of analysis used

    % disp(['OG index: ', num2str(index1)])

    % Get indices of duplicated elements
    dupeIdx = ismember(diffs,minDiff);
    dupeLoc = find(dupeIdx);

    % Logic Statements
    hasDuplicates = length(unique(rates)) < length(rates); % are there duplicates
    dup_is_min = any(dupeIdx == 1); % is the duplicate equal to the min difference

    if hasDuplicates && dup_is_min
        % if there are duplicates AND they equal the minimum difference, determine
        % min HR based off time to calculate
        [~ ,idx] = min(times(dupeLoc));
        index = dupeLoc(idx);
        HR = rates(index);
        % disp(['Time index: ', num2str(index)])
    else
        index = index1;
        HR = rates(index);
        % disp(['HR index: ', num2str(index)])
    end

end