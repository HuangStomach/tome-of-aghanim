
@Configuration
@EnableAsync
public class AsyncConfig {
    
    @Bean(name = "fileWriteExecutor")
    public Executor fileWriteExecutor() {
        ThreadPoolTaskExecutor executor = new ThreadPoolTaskExecutor();
        executor.setCorePoolSize(5);
        executor.setMaxPoolSize(10);
        executor.setQueueCapacity(100);
        executor.setThreadNamePrefix("FileWrite-");
        executor.initialize();
        return executor;
    }
}
