@echo off
title Run tests
jar cmf MainClass.txt submission.jar player99.class
javac -cp contest.jar player99.java Individual.java
jar -uf contest.jar Individual.class
echo Testing on Bent Cigar
java -jar testrun.jar -submission=player99 -evaluation=BentCigarFunction -seed=1
echo.
echo.
echo Testing on Schaffer's F7
java -jar testrun.jar -submission=player99 -evaluation=SchaffersEvaluation -seed=1
echo.
echo.
echo Testing on Katsuura
java -jar testrun.jar -submission=player99 -evaluation=KatsuuraEvaluation -seed=1
pause