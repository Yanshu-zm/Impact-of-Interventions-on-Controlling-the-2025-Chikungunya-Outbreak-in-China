function daysSinceBeginningOfYear = days_since_year(date_str)

specifiedDate = datetime(date_str); % Specify the date in the format 'yyyy-MM-dd'

% Get the beginning of the year for the specified date
beginningOfYear = datetime(year(specifiedDate), 1, 1);

% Calculate the number of days since the beginning of the year
daysSinceBeginningOfYear = days(specifiedDate - beginningOfYear);
end

