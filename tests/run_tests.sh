bash graphprot_test_misc.sh | grep --line-buffered -v "^Elapsed time" | sed -e 's/tmp-\S\{6\}/tmp-XXXXX/g' &> graphprot_test_misc.log

bash graphprot_test_classification_misc.sh | grep --line-buffered -v "^Elapsed time" | sed -e 's/tmp-\S\{6\}/tmp-XXXXX/g' &> graphprot_test_classification_misc.log

bash graphprot_test_regression_misc.sh | grep --line-buffered -v "^Elapsed time" | sed -e 's/tmp-\S\{6\}/tmp-XXXXX/g' &> graphprot_test_regression_misc.log

bash graphprot_test_regression_ls.sh | grep --line-buffered -v "^Elapsed time" | sed -e 's/tmp-\S\{6\}/tmp-XXXXX/g' &> graphprot_test_regression_ls.log

#bash graphprot_test_classification_ls.sh | grep --line-buffered -v "^Elapsed time" | sed -e 's/tmp-\S\{6\}/tmp-XXXXX/g' &> graphprot_test_classification_ls.log

gunzip -f *.gz