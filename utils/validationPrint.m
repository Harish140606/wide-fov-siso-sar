function ok = validationPrint(name, condition)
% VALIDATIONPRINT  Prints a "PASS"/"FAIL" line for a named check and
% returns the boolean so callers can aggregate an overall status.
% Only prints PASS if `condition` actually evaluated true.

if condition
    fprintf('  %-42s PASS\n', [name ':']);
else
    fprintf('  %-42s FAIL\n', [name ':']);
end
ok = logical(condition);
end
