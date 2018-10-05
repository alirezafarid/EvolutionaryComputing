javac -Xlint:unchecked  -Xlint:deprecation  -cp contest.jar player99.java individual.java CustomComparator.java
jar cmf MainClass.txt submission.jar player99.class individual.class CustomComparator.class 
java  -jar testrun.jar   -submission=player99  -evaluation=SchaffersEvaluation -seed=4
