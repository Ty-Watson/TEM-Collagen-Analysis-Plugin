
<h1 align="center" id="title">Collagen Fibril Analysis Plugin for ImageJ</h1>

<p align="center"><img src="https://socialify.git.ci/Ty-Watson/TEM-Collagen-Analysis-Plugin/image?font=Inter&amp;language=1&amp;name=1&amp;owner=1&amp;pattern=Signal&amp;stargazers=1&amp;theme=Light" alt="project-image"></p>

<p id="description">This ImageJ plugin is designed for the quantitative assessment of collagen fibril architecture from transmission electron microscopy (TEM) images inspired by the methodology described in the paper "A Fast Robust Method for Quantitative Assessment of Collagen Fibril Architecture from Transmission Electron Micrographs" by Bruno V. Rego Dar Weiss and Jay D. Humphrey.</p>

  
  
<h2>🧐 Features</h2>

Here're some of the project's best features:

*   Fibril Segmentation:
*   Quantitative Metrics
*   Interactive

<h2>🛠️ Installation Steps:</h2>

<p>1. Prereqs: Download and install ImageJ. Java Development Kit (JDK) - Version 8 or higher. Maven - A Java build tool you can install it from Maven's official site.</p>

<p>2. Clone the Repository First clone the plugin repository to your local machine. Open your terminal (or command prompt) and navigate to your desired directory and run:</p>

```
git clone https://github.com/Ty-Watson/TEM-Collagen-Analysis-Plugin.git
```

<p>3. Change to the directory where the plugin's source code is located:</p>

```
cd TEM-Collagen-Analysis-Plugin
```

<p>4. Use Maven to build the plugin into a JAR file. This will compile the Java files and package them into a plugin JAR file:</p>

```
mvn clean package
```

<p>5. Next you will copy the jar file to your ImageJ installation</p>

```
mvn -Dscijava.app.directory=/path/to/ImageJ.app/
```

<p>6. Run ImageJ and you should see the plugin "Collagen Analysis" under the Plugins Menu in ImageJ</p>

<p>7. Open an image with ImageJ then click the on the plugin in the menu to run</p>

<p>8. To update the plugin pull the latest changes from the repository and rebuild the JAR file:</p>

```
git pull
```

```
mvn clean package
```

  
  
<h2>💻 Built with</h2>

Technologies used in the project:

*   Java





