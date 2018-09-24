javac -cp contest.jar player99.java
jar cmf MainClass.txt submission.jar player99.class 
java  -jar testrun.jar   -submission=player99  -evaluation=BentCigarFunction -seed=1
