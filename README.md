
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

<h2>Manual copy</h2>
<p>1. Look into the target directory and find the generated jar file(tem-collagen-analysis-plugin-{version}.jar)</p>
<p>2. Drag and drop or copy the JAR file into the plugins folder in your ImageJ installation</p>
<p>3. Exit and reopen ImageJ and open the Plugins menu to verify the plugin "Collagen Analysis" is there</p>

```
mvn -Dscijava.app.directory=/path/to/ImageJ.app/
```

<p>6. Run ImageJ and you should see the plugin "Collagen Analysis" under the Plugins menu in ImageJ</p>

<p>7. Open an image with ImageJ then click the on the plugin in the menu to run</p>

<p>8. To update the plugin pull the latest changes from the repository and rebuild the JAR file:</p>

```
git pull
```

```
mvn clean package
```

<h2>Plugin Walkthrough</h2>
<p>1. Open an the image you want to process in ImageJ (file > open)</p>
<p>2. Find the Collagen Analysis plugin in the Plugins menu and click on it</p>
<p>3. Draw a line for the scale bar then click ok on the dialog box</p>
<p>4. Enter how long the line is in nanometers so the conversion factor can be set</p>
<p>5. Click in the image with the extended borders. If there are regions that need to be excluded from processing, then draw a polygon around the region you want to exclude then click a on you keyboard to exclude the region</p>
<p>6. When finished excluding regions in the image, press escape on your keyboard</p>
<p>7. Once the voronoi diagram is generated, right click on centroids that did not perform well to delete them. Left click areas where centroids need to be added. Once you are finished with editing the voronoi image, click ok in the dialog box</p>
<p>8. Once the image with the ellipses have been generated, left click on any ellipses that did not perform well to delete them. When you are finshed, press ok on the dialog box</p>
<p>9. Then, a file picker will prompt you to pick a directory to save the generated images and csv files. The folder name will be the name of the image you are processing</p>
<p>10. Look in the directory you specified to save the images to see the results</p>
  
<h2>💻 Built with</h2>

Technologies used in the project:

*   Java





