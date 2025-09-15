
package com.example.multidatasource;

import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import org.springframework.boot.autoconfigure.domain.EntityScan;
import org.springframework.data.jpa.repository.config.EnableJpaRepositories;

@SpringBootApplication
@EntityScan(basePackages = {
        "com.example.multidatasource.model.a",
        "com.example.multidatasource.model.b",
        "com.example.multidatasource.model.c"
})
@EnableJpaRepositories(basePackages = {
        "com.example.multidatasource.repository.a",
        "com.example.multidatasource.repository.b",
        "com.example.multidatasource.repository.c"
})
public class MultiDatasourceApplication {
    public static void main(String[] args) {
        SpringApplication.run(MultiDatasourceApplication.class, args);
    }
}
