function quantifyx(mitox, primaryB, layer2B, layer34B, terminalB, folderPath, xxshift, yyshift)
    %% quantify the size/ number/ circularity of mitochondria in each branches
[num_primary,size_primary_mean,circu_primary_mean, size_primary_all,circu_primary_all] = gen_score(mitox, primaryB, xxshift, yyshift);
fileID = fopen([folderPath,'mito_primary_branch_summary.txt'],'w');
fprintf(fileID,'number %d\n', num_primary);
fprintf(fileID,'mean size %f\n', size_primary_mean);
fprintf(fileID,'mean circularity %f\n', circu_primary_mean);
fprintf(fileID, ['ID ' 'size ' 'circularity\n']);
for i = 1:length(size_primary_all)
    fprintf(fileID, '%d %d %f\n', i,size_primary_all(i), circu_primary_all(i));
end
fclose(fileID);


[num_2,size_2_mean,circu_2_mean, size_2_all,circu_2_all] = gen_score(mitox, layer2B, xxshift, yyshift);
fileID = fopen([folderPath,'mito_secondary_branch_summary.txt'],'w');
fprintf(fileID,'number%d\n', num_2);
fprintf(fileID,'mean size %f\n', size_2_mean);
fprintf(fileID,'mean circularity %f\n', circu_2_mean);
fprintf(fileID, ['ID ' 'size ' 'circularity\n']);
for i = 1:length(size_2_all)
    fprintf(fileID, '%d %d %f\n', i,size_2_all(i), circu_2_all(i));
end
fclose(fileID);

[num_34,size_34_mean,circu_34_mean, size_34_all,circu_34_all] = gen_score(mitox, layer34B, xxshift, yyshift);
fileID = fopen([folderPath,'mito_fine_branch_summary.txt'],'w');
fprintf(fileID,'number%d\n', num_34);
fprintf(fileID,'mean size %f\n', size_34_mean);
fprintf(fileID,'mean circularity %f\n', circu_34_mean);
fprintf(fileID, ['ID ' 'size ' 'circularity\n']);
for i = 1:length(size_34_all)
    fprintf(fileID, '%d %d %f\n', i,size_34_all(i), circu_34_all(i));
end
fclose(fileID);

[num_terminal,size_terminal_mean,circu_terminal_mean, size_terminal_all,circu_terminal_all] = gen_score(mitox, terminalB, xxshift, yyshift);
fileID = fopen([folderPath,'mito_terminal_branch_summary.txt'],'w');
fprintf(fileID,'number%d\n', num_terminal);
fprintf(fileID,'mean size %f\n', size_terminal_mean);
fprintf(fileID,'mean circularity %f\n', circu_terminal_mean);
fprintf(fileID, ['ID ' 'size ' 'circularity\n']);
for i = 1:length(size_terminal_all)
    fprintf(fileID, '%d %d %f\n', i,size_terminal_all(i), circu_terminal_all(i));
end
fclose(fileID);

end