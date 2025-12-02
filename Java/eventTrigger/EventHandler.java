public class EventHandler {
	private boolean stopped = false;
	private int ttl = 128;
	private Object returnValue;

	public int ttl() {
		return this.ttl;
	}

	public EventHandler call() {
		this.ttl--;
		return this;
	}

	public boolean isStopped() {
		return this.stopped;
	}

	public void stopPropagation() {
		this.stopped = true;
	}

	public Object getReturnValue() {
		return returnValue;
	}

	public void setReturnValue(Object returnValue) {
		this.returnValue = returnValue;
	}
}
