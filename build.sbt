ThisBuild / version := "0.1.0-SNAPSHOT"

ThisBuild / scalaVersion := "2.13.10"

ThisBuild / scalacOptions := Seq("-deprecation", "-unchecked")

lazy val root = (project in file("."))
  .settings(
    name := "MultiViewLearningProject"
  )
libraryDependencies ++= Seq(
  "org.apache.spark" %% "spark-core" % "3.2.2",
  "org.apache.spark" %% "spark-sql" % "3.2.2",
  "org.scalanlp" %% "breeze" % "2.1.0",
  "org.scalanlp" %% "breeze-natives" % "2.1.0",
  "org.scalanlp" %% "breeze-viz" % "2.1.0",
  "org.scalatest" %% "scalatest" % "3.2.15" % "test"
)