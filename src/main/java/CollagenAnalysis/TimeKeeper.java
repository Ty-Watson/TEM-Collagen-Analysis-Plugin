package CollagenAnalysis;


public class TimeKeeper {
    private long startTime;
    private long endTime;

    // Start the timer
    public void start() {
        startTime = System.currentTimeMillis();
    }

    // Stop the timer and calculate the duration
    public void stop() {
        endTime = System.currentTimeMillis();
    }

    // Get the elapsed time in seconds
    public long getElapsedTimeSecs() {
        return (endTime - startTime) / 1000;
    }

    // Get the elapsed time in minutes and remaining seconds
    public String getElapsedTimeFormatted() {
        long duration = endTime - startTime;
        long minutes = (duration / 1000) / 60;
        long seconds = (duration / 1000) % 60;
        if (minutes > 0) {
            return String.format("%d minutes, %d seconds", minutes, seconds);
        } else {
            return String.format("%d seconds", seconds);
        }
    }
}
